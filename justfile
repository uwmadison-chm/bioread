# Run all checks
default: lint typecheck test

# Run the test suite; extra arguments go to pytest
test *args:
    uv run pytest {{args}}

# Run the tests against the lowest dependency versions we allow
test-lowest *args:
    uv run --isolated --resolution lowest-direct --all-extras pytest {{args}}

# Check linting and formatting
lint:
    uv run ruff check .
    uv run ruff format --check .

# Type check with ty
typecheck:
    uv run ty check

# Fix what ruff can fix automatically, and format
fix:
    uv run ruff check --fix .
    uv run ruff format .

# Build the sdist and wheel into dist/
build:
    uv build

# Update uv.lock after changing dependencies
lock:
    uv lock
