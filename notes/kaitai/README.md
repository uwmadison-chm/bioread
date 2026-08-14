# Describing AcqKnowledge headers in Kaitai Struct

An experiment: could bioread's header definitions and version-detection logic be
replaced by a [Kaitai Struct](https://kaitai.io/) spec, keeping streaming
support and getting something more robust than the current hand-rolled readers?

**Short answer: partly, and not for speed.** The header *layouts* map onto
Kaitai almost perfectly. The logic that *locates* those headers in a file does
not map at all, and would have to stay in Python. Kaitai parses headers roughly
50x slower than ctypes does.

Nothing here is wired into the bioread package. It is a spike, kept so the
next person to have this idea can start from the evidence.

## What's in here

| File | What it is |
| --- | --- |
| `gen_ksy.py` | Generates `acq.ksy` from the `_versioned_fields` tables in `bioread/headers.py` |
| `acq.ksy` | The generated spec (checked in so you can read it without a toolchain) |
| `compile_ksy.js` | Driver for the JavaScript build of the Kaitai compiler |
| `package.json` | Pins the compiler and a YAML parser |
| `validate.py` | Parses every file in `test/data` both ways and diffs the results |
| `benchmark.py` | Cost and bytes-touched comparison |

Reproducing it:

```sh
pip install kaitaistruct
npm install --prefix notes/kaitai

python notes/kaitai/gen_ksy.py > notes/kaitai/acq.ksy
node notes/kaitai/compile_ksy.js notes/kaitai/acq.ksy notes/kaitai
python notes/kaitai/validate.py
python notes/kaitai/benchmark.py
```

The spec is *generated* rather than hand-written on purpose. `headers.py`
already stores each layout as a table of `(name, ctype, minimum_revision)`,
which is very nearly a `.ksy` already, and generating means the two cannot
drift apart.

## What works, and better than expected

**Every graph header field, on every test file.** `validate.py` compares the
byte order, the file revision, and each graph header field the revision
includes:

```
files: 32   graph header fields compared: 1274
mismatches: 0
```

**The byte-order heuristic is expressible.** `GraphHeaderReader.bootstrap`
reads the version field both ways and keeps whichever gives the smaller number.
That works declaratively, so the spec orients itself:

```yaml
instances:
  version_le: {pos: 2, type: u4le}
  version_be: {pos: 2, type: u4be}
  is_le:         {value: version_le <= version_be}
  file_revision: {value: 'is_le ? version_le : version_be'}
```

Every type then resolves `endian: {switch-on: _root.is_le, ...}`.

**Version-gated fields map one-to-one.** A field's minimum revision becomes
`if: _root.file_revision >= 38`. Absent fields come back as `None`, which is
easier to work with than the current `hasattr` checks.

**Streaming is not a problem.** Positional `instances` are lazy — they parse
on attribute access, so a spec built from instances rather than `seq` never
touches the data block. Because the layout is lazy, a truncated or malformed
section also fails in isolation rather than taking the whole parse down, which
is roughly the behavior `Reader.read_errors` provides today.

**Python can drive it.** Top-level `params:` become constructor arguments, so
the offsets Python computes can be handed to Kaitai:

```python
parsed = Acq(data_offset, is_little_endian, KaitaiStream(f))
parsed.marker_header.l_markers   # seeks there, reads 8 bytes
```

## What doesn't work

**No `sum`, no `map`, no lambdas** in the expression language. Both fail to
compile:

```
error: don't know how to call method 'sum' of object type 'CalcArrayType(...)'
error: parsing expression 'items.map(x => x.len).sum' failed
```

So `data_length = Σ(point_count[i] × sample_size[i])` — the value that says
where markers begin in an uncompressed file — cannot live in the spec. The
top-level offset arithmetic stays in Python, and the "one spec describes the
whole file" idea does not survive.

**No scanning or backtracking.** `ChannelDTypeHeaderReader.scan_for_dtype_headers`
walks forward up to 4096 bytes looking for a plausible run of dtype headers.
Kaitai has no search primitive and no way to try-parse-and-retry, so this stays
in Python too.

**`encoding:` must be a compile-time literal.** bioread picks latin1 or utf-8
from the revision it just read, which cannot be expressed. `gen_ksy.py` emits
`size: N` for character arrays and leaves decoding to the caller. That is
arguably better layering, but it does mean a Python facade is still required —
you would not hand the generated class to users directly.

**The compiler is a JVM/JS program.** Not a runtime dependency (you check in
the generated parser), but it is a build step. The runtime, `kaitaistruct`, is
a small pure-Python package.

## Speed

Measured by `benchmark.py` on `physio-4.4.0.acq`:

```
Parsing one graph header:
  ctypes: for_revision + class + parse             48.2 us
  ctypes: build the struct class only              24.0 us
  ctypes: parse only, class reused                  0.5 us
  kaitai: construct + parse                        28.0 us
```

ctypes parsing is ~50x faster than Kaitai: `from_buffer_copy` is a memcpy,
while Kaitai does a `struct.unpack` per field. Kaitai would be a slowdown, not
a speedup.

The `build the struct class only` line is the interesting one. It used to
dominate, because `Header._struct_class` cached on the instance while
`for_revision()` returned a fresh instance per header — so every header read
rebuilt its ctypes class. That is fixed separately; it is not a reason to adopt
Kaitai.

One caveat on the generated spec: Kaitai reads a `repeat: expr` array one
element at a time, which is disastrous for the large opaque padding blobs in
these headers — an early version of this spec took 972 us per graph header and
issued 2280 read calls. `gen_ksy.py` emits `size: N` for byte arrays instead,
giving the 28 us above. If you hand-write a `.ksy`, do the same.

## Where this leaves things

Not worth a full migration. You would trade a Python file you can set a
breakpoint in for a YAML file, a JVM build step, and generated code you can't,
and you would give up most of the header-parsing speed — while still keeping
Python in charge of the offset arithmetic and the dtype scan, which are the
genuinely bodgy parts.

What is worth taking:

- **Publish `acq.ksy` as documentation.** It is generated and validated against
  the test corpus, so it stays honest, and it gives the format a machine-readable
  description plus readers in other languages for free. `validate.py` can run in
  CI as a check that the two descriptions still agree.
- **Steal the sizing idea.** Kaitai chains variable-length structures by reading
  a length and then a substream of that size. `reader.py` computes
  `foreign_offset` as `len(channel_headers) * ch_len`, assuming every channel
  header is the same length, while `HeaderReader.multi_headers` chains each
  header by its own `effective_len_bytes`. Those disagree if the lengths ever
  vary.
- **Note what became declarative.** The journal-magic rewind in
  `V2MarkerReader._read_marker_metadata` is a peeked tag and a conditional, and
  Kaitai's `--debug` mode records per-field byte offsets — which is what
  `acq_layout` reconstructs by hand.
