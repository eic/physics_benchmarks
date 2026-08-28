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
        float energy = energy_vec[i];
        float px = momx_vec[i];
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
    using ROOT::Experimental::RNTupleReader;
    
    auto ntuple = RNTupleReader::Open("events", rec_file);
    
    auto viewRecoNRG = ntuple->GetView<float>("ReconstructedChargedJets.energy");
    auto viewRecoMomZ = ntuple->GetView<float>("ReconstructedChargedJets.momentum.z");
    
    for (auto entryId : *ntuple) {
        auto energy_vec = viewRecoNRG(entryId);
        auto momz_vec = viewRecoMomZ(entryId);
        
        for (size_t i = 0; i < energy_vec.size(); i++) {
            float energy = energy_vec[i];
            float pz = momz_vec[i];
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

### RDataFrame Alternative

If the migration to native RNTupleReader is complex, consider using RDataFrame instead, which works with both TTree and RNTuple:

```cpp
// Works for both TTree and RNTuple
ROOT::RDataFrame df("events", rec_file);

auto df_filtered = df.Filter("ReconstructedChargedJets.energy.size() > 0")
                     .Define("jet_pt", "sqrt(ReconstructedChargedJets.momentum.x[0]*ReconstructedChargedJets.momentum.x[0] + "
                                       "ReconstructedChargedJets.momentum.y[0]*ReconstructedChargedJets.momentum.y[0])");
```

This approach requires minimal code changes but may have different performance characteristics.

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
