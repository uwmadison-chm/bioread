// Compile a .ksy spec to Python.
//
// The reference Kaitai compiler is a Scala program; this uses the JavaScript
// build of it, which only exposes a programmatic API, so we need this driver.
//
//   npm install --prefix notes/kaitai
//   node notes/kaitai/compile_ksy.js notes/kaitai/acq.ksy notes/kaitai
//
// The JVM command-line compiler from https://kaitai.io/#download does the same
// job if you would rather not have node in the loop.

const compiler = require('kaitai-struct-compiler');
const yaml = require('js-yaml');
const fs = require('fs');
const path = require('path');

const source = process.argv[2];
const outDir = process.argv[3] || '.';

if (!source) {
  console.error('usage: node compile_ksy.js <spec.ksy> [output-dir]');
  process.exit(2);
}

const spec = yaml.load(fs.readFileSync(source, 'utf8'));

compiler.compile('python', spec, null, false).then(
  files => {
    for (const name of Object.keys(files)) {
      const target = path.join(outDir, name);
      fs.writeFileSync(target, files[name]);
      console.log('wrote ' + target);
    }
  },
  err => {
    console.error('compile failed:', err.message || err);
    process.exit(1);
  }
);
