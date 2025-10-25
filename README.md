# ccmat

`ccmat` include subcrates, and all APIs are in crates `ccmat`. 

- `ccmat-core`: basic structure definition and utils, no other `ccmat-*` dependencies and all other crate depend on it.
- `ccmat-babel`: structure input parser to convert structure in between different format.
- `ccmat-kspace`: contain all reciprocal space operations. The `ckpath` is moved as independent crate. 

