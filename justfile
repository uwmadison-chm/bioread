# Run all checks
default: test

# Run the test suite; extra arguments go to pytest
test *args:
    uv run pytest {{args}}

# Build the sdist and wheel into dist/
build:
    uv build

# Update uv.lock after changing dependencies
lock:
    uv lock
