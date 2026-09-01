# RNTuple Migration Guide

This guide provides instructions for migrating physics benchmarks from TTree to RNTuple format.

## Overview

RNTuple is ROOT's next-generation columnar storage format, offering improved performance and better compression compared to TTree. This guide covers the necessary changes to migrate benchmark analysis code and Snakemake workflows.

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Simulation and Reconstruction Changes](#simulation-and-reconstruction-changes)
3. [Analysis Code Migration](#analysis-code-migration)
4. [Snakefile Updates](#snakefile-updates)
5. [Common Migration Patterns](#common-migration-patterns)
6. [Troubleshooting](#troubleshooting)

## Prerequisites

- ROOT 6.28 or later (RNTuple support)
- EICrecon with podio RNTuple backend support
- Updated EDM4eic/EDM4hep with RNTuple support

## Simulation and Reconstruction Changes

### EICrecon Output Configuration

To enable RNTuple output from eicrecon, use the podio backend configuration:

```bash
# Old TTree format (default)
eicrecon input.edm4hep.root -Ppodio:output_file=output.edm4eic.root

# New RNTuple format
eicrecon input.edm4hep.root -Ppodio:output_file=output.edm4eic.rnt.root -Ppodio:output_backend=rntuple
```

**Note:** Use `.rnt.root` extension for RNTuple files to distinguish them from TTree `.root` files.

### File Naming Convention

To distinguish RNTuple files from TTree files, use `.rnt.root` extension:

```
# TTree format
pythia8NCDIS_10x100_minQ2=1.edm4eic.root

# RNTuple format
pythia8NCDIS_10x100_minQ2=1.edm4eic.rnt.root
```

## Analysis Code Migration

### Overview of Changes

The main differences between TTree and RNTuple APIs:

| Feature | TTree API | RNTuple API |
|---------|-----------|-------------|
| Reader | `TTreeReader` | `ROOT::Experimental::RNTupleReader` |
| Field access | `TTreeReaderArray<T>` | `RNTupleReader::GetView<T>()` |
| Tree name | Required ("events") | Not used (single dataset per file) |
| Iteration | `reader.Next()` | Range-based for loop or manual iteration |

### Header Files

Replace TTree headers with RNTuple headers:

```cpp
// Old TTree includes
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>

// New RNTuple includes
#include <ROOT/RNTupleReader.hxx>
```

### Opening Files

#### TTree Approach (Old)

```cpp
TChain *mychain = new TChain("events");
mychain->Add(rec_file.c_str());
TTreeReader tree_reader(mychain);
```

#### RNTuple Approach (New)

```cpp
using ROOT::RNTupleReader;

auto ntuple = RNTupleReader::Open("events", rec_file);
if (!ntuple) {
  fmt::print(stderr, "ERROR: Failed to open RNTuple from file\n");
  return 1;
}
```

**Note:** In ROOT 6.40+, RNTuple is in the `ROOT::` namespace (not `ROOT::Experimental::`). Use `ROOT::RNTupleReader` for production code.

**Note:** RNTuple does not support chaining multiple files like TChain. If you need to process multiple files, you must:
1. Process them sequentially in a loop, or
2. Use RDataFrame which can handle multiple RNTuple files

### Error Handling

**Critical:** Unlike TTree which silently continues with missing branches, RNTuple throws exceptions for missing fields. Always wrap view creation in try-catch:

```cpp
try {
  auto viewRecoNRG = ntuple->GetView<float>("ReconstructedChargedJets.energy");
  auto viewRecoMomX = ntuple->GetView<float>("ReconstructedChargedJets.momentum.x");
  // ... more views
} catch (const std::exception& e) {
  fmt::print(stderr, "ERROR: Missing required field: {}\n", e.what());
  return 1;
}
```

This ensures your analysis fails fast if required collections are missing, rather than producing incorrect results silently.

### Reading Fields

#### TTree Approach (Old)

```cpp
TTreeReaderArray<float> recoNRG = {tree_reader, "ReconstructedChargedJets.energy"};
TTreeReaderArray<float> recoMomX = {tree_reader, "ReconstructedChargedJets.momentum.x"};

// In event loop
while (tree_reader.Next()) {
    for (int i = 0; i < recoNRG.GetSize(); i++) {
        float energy = recoNRG[i];
        float px = recoMomX[i];
        // ... process data
    }
}
```

#### RNTuple Approach (New)

```cpp
auto viewRecoNRG = ntuple->GetView<float>("ReconstructedChargedJets.energy");
auto viewRecoMomX = ntuple->GetView<float>("ReconstructedChargedJets.momentum.x");

for (auto entryId : *ntuple) {
    // For vector/array fields, you get a RVec-like object
    auto energy_vec = viewRecoNRG(entryId);
    auto momx_vec = viewRecoMomX(entryId);

    for (size_t i = 0; i < energy_vec.size(); i++) {
        float energy = energy_vec.at(i);  // Use .at() for ROOT compatibility
        float px = momx_vec.at(i);
        // ... process data
    }
}
```

**Alternative: Manual entry iteration**

```cpp
for (auto i = 0; i < ntuple->GetNEntries(); ++i) {
    auto energy_vec = viewRecoNRG(i);
    // ... process
}
```

### Field Type Mapping

When getting views, use the appropriate type:

```cpp
// Integer fields
auto viewType = ntuple->GetView<int>("ReconstructedChargedJets.type");
auto viewIndex = ntuple->GetView<int>("_ReconstructedChargedJets_constituents.index");

// Float fields
auto viewEnergy = ntuple->GetView<float>("ReconstructedChargedJets.energy");
auto viewMomX = ntuple->GetView<float>("ReconstructedChargedJets.momentum.x");

// Double fields (common in MCParticles)
auto viewMCMomX = ntuple->GetView<double>("MCParticles.momentum.x");

// Unsigned int fields
auto viewBegin = ntuple->GetView<unsigned int>("ReconstructedChargedJets.constituents_begin");
auto viewEnd = ntuple->GetView<unsigned int>("ReconstructedChargedJets.constituents_end");
```

**Note on Collection Indexing:**  
The exact syntax for indexing into collections returned by `view(entryId)` may vary by ROOT version. If you encounter errors like "subscripted value is not an array", try these alternatives:
- Use `.at(i)` instead of `[i]`
- Use range-based for loops: `for (const auto& val : view(entryId))`
- Check your ROOT version's RNTuple documentation for the exact API

### Complete Migration Example

Here's a complete before/after example:

#### Before (TTree)

```cpp
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>

int analyze(const std::string& rec_file) {
    TChain *mychain = new TChain("events");
    mychain->Add(rec_file.c_str());
    TTreeReader tree_reader(mychain);
    
    TTreeReaderArray<float> recoNRG = {tree_reader, "ReconstructedChargedJets.energy"};
    TTreeReaderArray<float> recoMomZ = {tree_reader, "ReconstructedChargedJets.momentum.z"};
    
    while (tree_reader.Next()) {
        for (int i = 0; i < recoNRG.GetSize(); i++) {
            float energy = recoNRG[i];
            float pz = recoMomZ[i];
            // Process jet...
        }
    }
    return 0;
}
```

#### After (RNTuple)

```cpp
#include <ROOT/RNTupleReader.hxx>

int analyze(const std::string& rec_file) {
    using ROOT::RNTupleReader;  // ROOT 6.40+: use ROOT:: not ROOT::Experimental::
    
    auto ntuple = RNTupleReader::Open("events", rec_file);
    
    auto viewRecoNRG = ntuple->GetView<float>("ReconstructedChargedJets.energy");
    auto viewRecoMomZ = ntuple->GetView<float>("ReconstructedChargedJets.momentum.z");
    
    for (auto entryId : *ntuple) {
        auto energy_vec = viewRecoNRG(entryId);
        auto momz_vec = viewRecoMomZ(entryId);

        for (size_t i = 0; i < energy_vec.size(); i++) {
            float energy = energy_vec.at(i);  // Use .at() for ROOT compatibility
            float pz = momz_vec.at(i);
            // Process jet...
        }
    }
    return 0;
}
```

## Snakefile Updates

### Update Reconstruction Rule

Modify the reconstruction rule to use RNTuple output:

```python
# Before
rule my_reco_eicrecon:
    input:
        "sim_output/{DETECTOR_CONFIG}/input.edm4hep.root",
    output:
        "sim_output/{DETECTOR_CONFIG}/output.edm4eic.root",
    shell:
        """
        DETECTOR_CONFIG={wildcards.DETECTOR_CONFIG} eicrecon {input} -Ppodio:output_file={output}
        """

# After
rule my_reco_eicrecon:
    input:
        "sim_output/{DETECTOR_CONFIG}/input.edm4hep.root",
    output:
        "sim_output/{DETECTOR_CONFIG}/output.edm4eic.rnt.root",
    shell:
        """
        DETECTOR_CONFIG={wildcards.DETECTOR_CONFIG} eicrecon {input} \\
            -Ppodio:output_file={output} \\
            -Ppodio:output_backend=rntuple
        """
```

### Update File References

Update all file references throughout the Snakefile:

```python
# Before
data="sim_output/{DETECTOR_CONFIG}/pythia8NCDIS_10x100.edm4eic.root",

# After
data="sim_output/{DETECTOR_CONFIG}/pythia8NCDIS_10x100.edm4eic.rnt.root",
```

## Common Migration Patterns

### Pattern 1: Simple Field Access

```cpp
// Old
TTreeReaderArray<float> field = {tree_reader, "Collection.field"};
while (tree_reader.Next()) {
    float value = field[index];
}

// New
auto viewField = ntuple->GetView<float>("Collection.field");
for (auto entryId : *ntuple) {
    auto field_vec = viewField(entryId);
    float value = field_vec.at(index);  // Use .at() for safety and ROOT compatibility
}
```

### Pattern 2: Checking Array Size

```cpp
// Old
int size = recoNRG.GetSize();

// New
auto energy_vec = viewRecoNRG(entryId);
size_t size = energy_vec.size();
```

### Pattern 3: Conditional Field Access (Version-Dependent)

```cpp
// Old
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8,9,0)
  TTreeReaderArray<float> recoArea = {tree_reader, "ReconstructedChargedJets.area"};
#endif

// New
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8,9,0)
  auto viewRecoArea = ntuple->GetView<float>("ReconstructedChargedJets.area");
#endif
```

### Pattern 4: Association Tables

```cpp
// Old
TTreeReaderArray<unsigned int> recoCstsBegin = {tree_reader, "ReconstructedChargedJets.constituents_begin"};
TTreeReaderArray<unsigned int> recoCstsEnd = {tree_reader, "ReconstructedChargedJets.constituents_end"};
TTreeReaderArray<int> recoCstIndex = {tree_reader, "_ReconstructedChargedJets_constituents.index"};

while (tree_reader.Next()) {
    for (int i = 0; i < recoType.GetSize(); i++) {
        unsigned int begin = recoCstsBegin[i];
        unsigned int end = recoCstsEnd[i];
        for (unsigned int j = begin; j < end; j++) {
            int idx = recoCstIndex[j];
            // Process constituent at idx
        }
    }
}

// New
auto viewRecoCstsBegin = ntuple->GetView<unsigned int>("ReconstructedChargedJets.constituents_begin");
auto viewRecoCstsEnd = ntuple->GetView<unsigned int>("ReconstructedChargedJets.constituents_end");
auto viewRecoCstIndex = ntuple->GetView<int>("_ReconstructedChargedJets_constituents.index");

for (auto entryId : *ntuple) {
    auto cstsBegin = viewRecoCstsBegin(entryId);
    auto cstsEnd = viewRecoCstsEnd(entryId);
    auto cstIndex = viewRecoCstIndex(entryId);
    
    for (size_t i = 0; i < cstsBegin.size(); i++) {
        unsigned int begin = cstsBegin[i];
        unsigned int end = cstsEnd[i];
        for (unsigned int j = begin; j < end; j++) {
            int idx = cstIndex[j];
            // Process constituent at idx
        }
    }
}
```

## Troubleshooting

### "Cannot open RNTuple" Error

**Cause:** Trying to open a TTree file with RNTupleReader or vice versa.

**Solution:** Ensure the file was created with `-Ppodio:output_backend=rntuple` and verify the format:

```bash
root -l -q 'file.root' -e 'gDirectory->ls()'
```

Look for `RNTuple` instead of `TTree` in the output.

### Type Mismatch Errors

**Cause:** Using wrong type in `GetView<T>()`.

**Solution:** Check the field type in the RNTuple:

```cpp
ntuple->GetDescriptor().PrintInfo();
```

Common types:
- `int` for PDG codes, type fields
- `float` for most EDM4eic fields (energy, momentum)
- `double` for MCParticles
- `unsigned int` for indices and counts

### Performance Issues

**Cause:** Creating views inside the event loop.

**Solution:** Always create views once before the loop:

```cpp
// BAD - creates view every iteration
for (auto entryId : *ntuple) {
    auto view = ntuple->GetView<float>("field");  // DON'T DO THIS
}

// GOOD - creates view once
auto view = ntuple->GetView<float>("field");
for (auto entryId : *ntuple) {
    auto data = view(entryId);
}
```

### Missing Fields

**Cause:** Field name changed or doesn't exist in RNTuple.

**Solution:** List all available fields:

```cpp
ntuple->GetDescriptor().PrintInfo();
// or
for (const auto& field : ntuple->GetDescriptor().GetFieldRange()) {
    std::cout << field.GetFieldName() << std::endl;
}
```

### podio DataFrame API (Recommended)

**The recommended approach** for EDM4hep/EDM4eic analysis is using podio's DataFrame API, which provides format-agnostic access and handles podio collections natively:

```cpp
#include <podio/DataSource.h>
#include <ROOT/RDataFrame.hxx>

int analyze(const std::string& rec_file) {
    // Automatically detects TTree vs RNTuple format
    auto df = podio::CreateDataFrame(rec_file);
    
    // Use standard RDataFrame operations
    // Collections are accessible with dot notation: "Collection.field"
    
    // Example: For complex per-event logic, use Foreach
    df.Foreach([&histogram1, &histogram2](
        ROOT::VecOps::RVec<float> jet_energy,
        ROOT::VecOps::RVec<float> jet_px,
        ROOT::VecOps::RVec<float> jet_py
    ) {
        // Your analysis logic with full control
        for (size_t i = 0; i < jet_energy.size(); i++) {
            float energy = jet_energy.at(i);
            float pt = sqrt(jet_px.at(i)*jet_px.at(i) + jet_py.at(i)*jet_py.at(i));
            histogram1->Fill(energy);
            histogram2->Fill(pt);
        }
    }, {
        "ReconstructedChargedJets.energy",
        "ReconstructedChargedJets.momentum.x",
        "ReconstructedChargedJets.momentum.y"
    });
    
    return 0;
}
```

**Advantages:**
- ✅ Format-agnostic: Works with `.root` (TTree) and `.rnt.root` (RNTuple) transparently
- ✅ Uses podio's official API for EDM4hep/EDM4eic data
- ✅ No manual view creation or error handling needed
- ✅ Cleaner initialization code
- ✅ Can be parallelized with `ROOT::EnableImplicitMT()` (ensure side effects like manual `TH1::Fill` are made thread-safe, e.g. via `ForeachSlot` + per-slot histograms)
- ✅ Perfect for complex analysis with nested constituent loops

**When to use `.Foreach()` vs pure declarative style:**
- Use `.Foreach()` for complex per-event logic (constituent loops, jet matching, multi-histogram filling)
- Use `.Define()` and `.Filter()` for simple transformations and cuts
- See `benchmarks/Jets-HF/jets/analysis/jets.cxx` for a complete real-world example

**Migration from RNTupleReader to podio::CreateDataFrame:**

```cpp
// Before (RNTupleReader)
#include <ROOT/RNTupleReader.hxx>

auto ntuple = RNTupleReader::Open("events", rec_file);
auto viewEnergy = ntuple->GetView<float>("ReconstructedChargedJets.energy");
auto viewMomX = ntuple->GetView<float>("ReconstructedChargedJets.momentum.x");

for (auto entryId : *ntuple) {
    auto energy = viewEnergy(entryId);
    auto momx = viewMomX(entryId);
    // ... process
}

// After (podio::CreateDataFrame)
#include <podio/DataSource.h>
#include <ROOT/RDataFrame.hxx>

auto df = podio::CreateDataFrame(rec_file);

df.Foreach([&histograms](
    ROOT::VecOps::RVec<float> energy,
    ROOT::VecOps::RVec<float> momx
) {
    // Same processing logic
}, {
    "ReconstructedChargedJets.energy",
    "ReconstructedChargedJets.momentum.x"
});
```

### Plain RDataFrame Alternative

If podio DataFrame is not available, you can use plain RDataFrame which also works with both TTree and RNTuple:

```cpp
// Works for both TTree and RNTuple (but less podio-aware)
ROOT::RDataFrame df("events", rec_file);

auto df_filtered = df.Filter("ReconstructedChargedJets.energy.size() > 0")
                     .Define("jet_pt", "sqrt(ReconstructedChargedJets.momentum.x[0]*ReconstructedChargedJets.momentum.x[0] + "
                                       "ReconstructedChargedJets.momentum.y[0]*ReconstructedChargedJets.momentum.y[0])");
```

This approach requires minimal code changes but may have different performance characteristics.

## Migration Considerations and Limitations

### Multiple File Processing (TChain Replacement)

**Limitation:** RNTuple does not have a direct equivalent to TChain for processing multiple files transparently.

**Solutions:**

1. **Sequential Processing** (simplest for benchmarks):
```cpp
for (const auto& filename : input_files) {
    auto ntuple = RNTupleReader::Open("events", filename);
    // Process each file
}
```

2. **RDataFrame with Multiple Files** (recommended for analysis):
```cpp
ROOT::RDataFrame df("events", {"file1.rnt.root", "file2.rnt.root", "file3.rnt.root"});
// Declarative analysis works across all files
```

3. **Format-Agnostic Wrapper** (for mixed TTree/RNTuple workflows):
Create a wrapper class that detects file format and uses TChain for TTree or sequential RNTupleReader for RNTuple.

### Backwards Compatibility Strategy

**Challenge:** Analyzing old TTree data alongside new RNTuple data requires either:
- Two versions of analysis code, or
- Format-agnostic code using RDataFrame

**Recommendations:**

1. **For Controlled Environments (Benchmarks)**:
   - Each campaign uses one format consistently
   - Cross-campaign comparisons can regenerate old data in RNTuple format if needed
   - This migration approach is suitable

2. **For General Analysis (Ongoing Studies)**:
   - Use RDataFrame which handles both formats transparently
   - Or maintain a thin abstraction layer that dispatches to TTreeReader or RNTupleReader based on file format

3. **Hybrid Approach**:
   - Store format detection logic in a utility function
   - Branch analysis code based on detected format
   - Example:
```cpp
bool isRNTuple(const std::string& filename) {
    return filename.find(".rnt.root") != std::string::npos ||
           filename.find(".rntuple.root") != std::string::npos;
}
```

### Python Support

**Status:** Both uproot (5.x) and ROOT's Python bindings (PyROOT) support RNTuple starting with ROOT 6.28+.

## Python Script Migration Patterns

Python analysis scripts can be migrated to work transparently with both TTree and RNTuple formats. The recommended approach is using **uproot 5.x**, which provides automatic format detection and a unified API for both formats.

### Overview of Python Migration Approaches

| Approach | Pros | Cons | Best For |
|----------|------|------|----------|
| **uproot 5.x** (Recommended) | ✅ Automatic format detection<br>✅ Pure Python, no ROOT install needed<br>✅ Excellent performance<br>✅ Clean, Pythonic API | ⚠️ Requires uproot >= 4.0 | Most Python analysis scripts |
| **PyROOT with podio** | ✅ Direct access to podio API<br>✅ Uses same pattern as C++ | ⚠️ Requires ROOT install<br>⚠️ Less Pythonic | Scripts that need C++ ROOT features |

### Format-Agnostic Python with uproot (Recommended)

**Key Insight:** uproot 5.x automatically handles both TTree and RNTuple formats with the same syntax. Podio creates RNTuple files with the same `events` tree/RNTuple name and branch structure as TTree files, so existing uproot code often works without modification.

#### Basic Pattern

```python
import uproot

# This works transparently with both:
# - .edm4eic.root (TTree format)
# - .edm4eic.rnt.root (RNTuple format)
keys = uproot.concatenate(rec_file + ':events/' + 'CollectionName')
data_Q2 = keys['CollectionName.Q2']
data_x = keys['CollectionName.x']
```

**Why it works:**
- uproot 5.x detects whether `events` is a TTree or RNTuple automatically
- Both formats use the same branch/field naming conventions from podio
- No code changes needed if branch names are consistent

#### Adding Format-Agnostic Comments

To document the format-agnostic capability and help future maintainers:

```python
# Format-agnostic data loading: uproot 5.x automatically handles both TTree and RNTuple formats
# Works with:
#   - .edm4eic.root files (TTree format)
#   - .edm4eic.rnt.root files (RNTuple format)
# Both formats created by podio use the same 'events' tree/RNTuple name and branch structure
keys = uproot.concatenate(rec_file + ':events/' + 'InclusiveKinematicsTruth')
Truth = [keys['InclusiveKinematicsTruth.Q2'], keys['InclusiveKinematicsTruth.x']]
```

#### Complete Real-World Example

See `benchmarks/Inclusive/dis/analysis/kinematics_correlations.py` for a complete working example:

```python
#!/usr/bin/env python
import uproot as ur
import awkward as ak
import numpy as np

# No format detection needed - uproot handles it automatically
rec_file = "data.edm4eic.root"  # or data.edm4eic.rnt.root

# Load data - works with both TTree and RNTuple
keys = ur.concatenate(rec_file + ':events/' + 'InclusiveKinematicsTruth')
Truth = [keys['InclusiveKinematicsTruth.Q2'], keys['InclusiveKinematicsTruth.x']]

keys = ur.concatenate(rec_file + ':events/' + 'InclusiveKinematicsElectron')
Electron = [keys['InclusiveKinematicsElectron.Q2'], keys['InclusiveKinematicsElectron.x']]

# Process data with awkward arrays (same code for both formats)
Q2values_T = Truth[0]
Xvalues_T = Truth[1]
T_Q2s = np.array(ak.flatten(Q2values_T))
T_Xs = np.array(ak.flatten(Xvalues_T))
```

### Advanced: Explicit Format Detection (Optional)

If you need to handle formats differently or provide informative logging:

```python
import uproot

def detect_format(filename):
    """Detect if file contains TTree or RNTuple."""
    with uproot.open(filename) as file:
        # Check what 'events' is
        events = file['events']
        if 'TTree' in events.classname:
            return 'TTree'
        elif 'RNTuple' in events.classname or hasattr(events, 'file'):
            return 'RNTuple'
    return 'Unknown'

# Example usage
rec_file = "data.edm4eic.root"
format_type = detect_format(rec_file)
print(f"Detected format: {format_type}")

# Same code works regardless of format
keys = uproot.concatenate(rec_file + ':events/' + 'Collection')
```

### Alternative: PyROOT with podio.CreateDataFrame

For scripts that need ROOT features or want to mirror C++ patterns:

```python
import ROOT
from podio import CreateDataFrame

# Load file (format-agnostic, just like C++)
rec_file = "data.edm4eic.root"  # or .rnt.root
df = CreateDataFrame(rec_file)

# Use RDataFrame operations
# Define columns
df = df.Define("Q2_truth", "InclusiveKinematicsTruth.Q2")
df = df.Define("x_truth", "InclusiveKinematicsTruth.x")

# Convert to numpy for plotting
Q2_array = df.AsNumpy(["Q2_truth"])["Q2_truth"]
x_array = df.AsNumpy(["x_truth"])["x_truth"]
```

**Pros:**
- Mirrors C++ approach exactly
- Access to full ROOT ecosystem
- Type-safe column operations

**Cons:**
- Requires ROOT installation
- Less Pythonic than uproot
- Harder to use with awkward arrays

### Best Practices for Format-Agnostic Python

1. **Use uproot 5.x or later**
   ```bash
   pip install 'uproot>=5.0'
   ```

2. **Don't hardcode format assumptions**
   ```python
   # ❌ Bad - assumes TTree
   tree = file['events']
   assert isinstance(tree, uproot.TTree)
   
   # ✅ Good - works with both
   events = file['events']
   data = events.arrays(['branch1', 'branch2'])
   ```

3. **Rely on consistent naming**
   - Podio ensures both TTree and RNTuple files have the same branch/field names
   - Use the same collection and field access patterns for both formats

4. **Add format-agnostic comments**
   - Document that code works with both formats
   - List the file extensions supported
   - Note any assumptions about branch structure

5. **Test with TTree first, RNTuple later**
   - Existing TTree files are the compatibility baseline
   - RNTuple files should work with the same code
   - If RNTuple doesn't work, that's a bug in the migration (not expected with uproot 5.x)

6. **Handle legacy compatibility explicitly**
   ```python
   # Handle renamed collections gracefully
   try:
       keys = ur.concatenate(rec_file + ':events/' + 'NewCollectionName')
   except ur.KeyInFileError:
       # Fallback for older files
       keys = ur.concatenate(rec_file + ':events/' + 'OldCollectionName')
   ```

### Migration Checklist for Python Scripts

- [ ] Verify uproot version is >= 5.0 (`import uproot; print(uproot.__version__)`)
- [ ] Add format-agnostic comment at data loading section
- [ ] Test script runs without errors (syntax, imports, basic logic)
- [ ] Verify awkward array operations work (no format-specific assumptions)
- [ ] Document any collection name fallbacks for legacy compatibility
- [ ] Note in comments that script works with `.edm4eic.root` and `.edm4eic.rnt.root`

### Troubleshooting Python Migration

#### "KeyInFileError: not found in file"

**Cause:** Branch/field name doesn't exist or collection name changed.

**Solution:**
```python
# List available collections
with uproot.open(rec_file) as file:
    print(file['events'].keys())  # Works for both TTree and RNTuple
```

#### "Different array lengths" or shape mismatches

**Cause:** Usually not format-related, but rather different event content.

**Solution:** Verify the input file was created with correct simulation/reconstruction parameters.

#### Performance differences between TTree and RNTuple

**Expected:** RNTuple may be faster for columnar access, especially for large files.

**Action:** No code changes needed. Document observed performance if significant.

### Summary

**For most Python analysis scripts:**
1. Ensure uproot >= 5.0 is installed
2. Add a comment documenting format-agnostic capability
3. Existing code using `uproot.concatenate(file + ':events/' + 'Collection')` should work unchanged
4. Test with existing TTree files to verify no regressions

**The key insight:** uproot 5.x + podio's consistent naming makes most Python migrations trivial - often just adding documentation rather than changing code.

For a complete working example, see the migration of `benchmarks/Inclusive/dis/analysis/kinematics_correlations.py`.

**Usage:**
```python
import ROOT

# Open RNTuple file
ntuple = ROOT.RNTupleReader.Open("events", "output.rnt.root")

# Access data
for entry in ntuple:
    energy_view = ntuple.GetView[float]("ReconstructedChargedJets.energy")
    # Process data
```

**Note:** Python analysis may have different ergonomics than C++. For production Python analysis, consider:
- Using RDataFrame Python interface for declarative analysis
- uproot library (check RNTuple support status in your uproot version)

## References

- [ROOT RNTuple Documentation](https://root.cern/doc/master/md_tree_ntuple_v7_doc_README.html)
- [RNTuple Tutorial](https://root.cern/doc/master/ntpl001__staff_8C.html)
- [EDM4eic Documentation](https://github.com/eic/EDM4eic)
- [Podio Documentation](https://github.com/AIDASoft/podio)

## Getting Help

If you encounter issues during migration:
1. Check this guide's troubleshooting section
2. Verify ROOT version supports RNTuple (>= 6.28)
3. Confirm podio backend is compiled with RNTuple support
4. Ask in the EIC Software Group Slack channel
5. Open an issue in the physics_benchmarks repository

## Migration Checklist

When migrating a benchmark:

- [ ] Update eicrecon command with `-Ppodio:output_backend=rntuple`
- [ ] Update file extensions to `.edm4eic.rnt.root`
- [ ] Replace `TChain`/`TTreeReader` includes with `RNTupleReader`
- [ ] Convert file opening to `RNTupleReader::Open()` with null check
- [ ] Add try-catch around all `GetView<T>()` calls for error handling
- [ ] Replace all `TTreeReaderArray` with `GetView<T>()`
- [ ] Update event loop from `while (reader.Next())` to `for (auto entryId : *ntuple)`
- [ ] Update array access from `field[i]` to `view(entryId).at(i)`
- [ ] Move view creation outside event loop for performance
- [ ] Test compilation
- [ ] Verify output produces expected results
- [ ] Update documentation/comments mentioning TTree
