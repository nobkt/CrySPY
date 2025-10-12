# Fix for Infinite Retry Loop When Generating Crystal Structures

## Problem

When PyXtal cannot generate a crystal structure for certain space groups, it outputs "Cannot generate crystal after max attempts." repeatedly, causing very long wait times or even infinite loops.

Example:
```
[2025-10-12 19:49:33,191][cryspy_mol_crystal][INFO] 空間群 143 で結晶生成を試行中...
Cannot generate crystal after max attempts.
Cannot generate crystal after max attempts.
Cannot generate crystal after max attempts.
...
```

## Solution

We've added retry limits and smart failure tracking to prevent infinite loops:

### 1. Maximum Attempt Limits

All structure generation functions now have a `maxcnt` parameter (default: 100):
- `gen_struc()` 
- `gen_struc_mol()`
- `gen_struc_mol_break_sym()`

### 2. Space Group Failure Tracking

The code now:
- Tracks how many times each space group fails
- Skips space groups that fail more than 10 times
- Provides clear error messages when max attempts reached

### 3. PyXtal Failure Detection

The molecular crystal script now:
- Captures PyXtal's stdout/stderr output
- Detects "Cannot generate crystal after max attempts" messages
- Skips problematic space groups and volume factors early

## Usage

### Using gen_pyxtal functions

```python
from cryspy.RS.gen_struc_RS.gen_pyxtal import gen_struc

try:
    structures = gen_struc(
        nstruc=10,
        atype=('Si',),
        nat=(8,),
        mindist=((2.0,),),
        spgnum=[1, 2, 14, 15],
        maxcnt=200,  # Increase if needed (default: 100)
    )
except RuntimeError as e:
    print(f"Failed to generate structures: {e}")
    # The error message will suggest adjusting parameters
```

### Using cryspy-mol-crystal script

The script automatically handles failures with sensible defaults:

```bash
cryspy-mol-crystal molecule.xyz -o output.cif --optimize-density
```

If a particular space group keeps failing, the script will:
1. Try multiple volume factors
2. Skip the space group after too many failures
3. Move on to the next space group
4. Report which space groups succeeded

## Adjusting Parameters

If you're having trouble generating structures, try:

1. **Increase max attempts**: Add `maxcnt=200` (or higher) parameter
2. **Adjust vol_factor**: Try different values like 1.0, 1.2, 1.5
3. **Relax mindist**: Use larger minimum distances
4. **Try different space groups**: Some space groups are easier than others
5. **Reduce structure complexity**: Try fewer atoms/molecules first

## Error Messages

When maximum attempts are reached, you'll see:
```
[ERROR] Reached maximum attempts (100) without generating a valid structure.
[ERROR] Space group failure counts: {143: 11, 144: 8, 145: 7, ...}
RuntimeError: Cannot generate structure after 100 attempts. 
Consider adjusting vol_factor, mindist, or space group selection.
```

This tells you:
- Which space groups failed and how many times
- Suggestions for what to adjust

## Technical Details

### Changes in `gen_pyxtal.py`

```python
# Before: Could loop forever
while len(init_struc_data) < nstruc:
    # ... try to generate structure ...
    # If fails, continue (infinite loop possible)

# After: Limited attempts with tracking
while len(init_struc_data) < nstruc:
    cnt = 0
    spg_fail_count = {}
    
    while cnt < maxcnt:
        cnt += 1
        # Skip space groups that failed too many times
        if spg_fail_count.get(spg, 0) > 10:
            continue
        # ... try to generate structure ...
        # Track failures per space group
        spg_fail_count[spg] += 1
        
    if cnt >= maxcnt:
        raise RuntimeError(...)  # Clear error message
```

### Changes in `cryspy_mol_crystal.py`

```python
# Capture PyXtal output
from contextlib import redirect_stdout, redirect_stderr
from io import StringIO

f = StringIO()
with redirect_stdout(f):
    with redirect_stderr(f):
        crystal.from_random(...)

output = f.getvalue()
if "Cannot generate crystal after max attempts" in output:
    # Track and skip this space group/volume factor
    failed_count += 1
```

## Backwards Compatibility

All changes are backwards compatible:
- Default parameter values maintain existing behavior
- Existing code will work without modification
- New parameters are optional with sensible defaults

## Testing

The changes have been validated with:
- Syntax checks
- Import tests  
- Parameter signature tests
- Smoke tests for all modified functions

For more information, see the code comments in:
- `src/cryspy/RS/gen_struc_RS/gen_pyxtal.py`
- `src/cryspy/scripts/cryspy_mol_crystal.py`
