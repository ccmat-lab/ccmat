use pyo3::prelude::*;

#[pymodule]
#[pyo3(name = "ccmat")]
fn ccmatpy(m: &Bound<'_, PyModule>) -> PyResult<()> {
    Ok(())
}
