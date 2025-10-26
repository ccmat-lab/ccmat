/// notes:
/// Operation safety is guranteed by the type.
/// Use Angstrom as the major internal and default API unit to be consistent with xx/xx.
/// Internally use FracCoord to represent the position to make sure fractional
/// coords don’t change if lattice changes shape or scale.
///
/// - Use 'lattice'

// TODO: naming convention for vars, check IUCr or cif specification
// Give a table to compare in between different popular tools.
//
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Angstrom(pub f64);

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bohr(pub f64);

// FracCoord [0.0, 1.0) and only used internally for site position.
// TODO: internally I need to check the Frac is valid in between 0.0~1.0
#[derive(Debug, Clone, Copy, PartialEq)]
struct FracCoord(f64);

/// Lattice
/// inner data structure of the struct are private.
#[derive(Debug)]
pub struct Lattice {
    a: [Angstrom; 3],
    b: [Angstrom; 3],
    c: [Angstrom; 3],
}

impl Lattice {
    // TODO: how to use type system to validate the row/column definition or unit?
    pub fn new(a: [Angstrom; 3], b: [Angstrom; 3], c: [Angstrom; 3]) -> Self {
        Lattice { a, b, c }
    }
}

// TODO: add lattice_bohr!()
// TODO: impl Display to Lattice for pretty print.

#[macro_export]
macro_rules! __vec3_angstrom {
    ([$x:expr, $y:expr, $z:expr]) => {
        [Angstrom($x), Angstrom($y), Angstrom($z)]
    };
    (($x:expr, $y:expr, $z:expr)) => {
        [
            $crate::Angstrom($x),
            $crate::Angstrom($y),
            $crate::Angstrom($z),
        ]
    };
}

/// Create a [`Lattice`] using vectors expressed in **Ångström** units.
///
/// This macro constructs a [`Lattice`] from three lattice vectors (`a`, `b`, and `c`)
/// expressed as tuples or arrays of three floating-point numbers.
///
/// - Each component is converted to [`Angstrom`] automatically.
/// - Both `(x, y, z)` and `[x, y, z]` tuple/array syntax are supported.
/// - Trailing commas are optional.
///
/// It supports both **named** and **positional** forms:
///
/// - **Named form** (explicit `a=`, `b=`, `c=`):
///   ```
///   use commat::lattice_angstrom;
///
///   let latt = lattice_angstrom!(
///       a = (1.0, 0.0, 0.0),
///       b = (0.0, 1.0, 0.0),
///       c = (0.0, 0.0, 1.0),
///   );
///   ```
///
/// - **Positional form** (omit names, ordered as `a`, `b`, `c`):
///   ```
///   use commat::{lattice_angstrom, Lattice, Angstrom};
///
///   let latt = lattice_angstrom!(
///       (1.0, 0.0, 0.0),
///       (0.0, 1.0, 0.0),
///       (0.0, 0.0, 1.0),
///   );
///   ```
///
/// # Errors
/// - None at compile time; this macro expands directly to constructor calls.
///
/// # Example
/// ```
/// use commat::{lattice_angstrom, Lattice, Angstrom};
///
/// let latt = lattice_angstrom!(
///     a = [2.5, 0.0, 0.0],
///     b = [0.0, 2.5, 0.0],
///     c = [0.0, 0.0, 2.5],
/// );
/// println!("{:?}", latt);
/// ```
#[macro_export]
macro_rules! lattice_angstrom {
    (
        a = $a:tt,
        b = $b:tt,
        c = $c:tt $(,)?
    ) => {{
        let lattice = $crate::Lattice::new(
            $crate::__vec3_angstrom!($a),
            $crate::__vec3_angstrom!($b),
            $crate::__vec3_angstrom!($c),
        );
        lattice
    }};
    (
        $a:tt,
        $b:tt,
        $c:tt $(,)?
    ) => {{
        let lattice = Lattice::new(
            $crate::__vec3_angstrom!($a),
            $crate::__vec3_angstrom!($b),
            $crate::__vec3_angstrom!($c),
        );
        lattice
    }};
}

#[derive(Debug)]
pub struct CrystalValidateError {
    message: String,
}

impl std::fmt::Display for CrystalValidateError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl std::error::Error for CrystalValidateError {}

pub struct LatticeSet;
pub struct LatticeNotSet;
pub struct AtomsSet;
pub struct AtomsNotSet;

/// use builder pattern so the validation is runtime
/// To make it a compile time check, I can use the proc macro (`build`/`build_unchcek` API).
///
/// # Example
/// ```
/// use commat::*;
///
/// let lattice = lattice_angstrom![
///     a = (1.0, 0.0, 0.0),
///     b = (0.0, 1.0, 0.0),
///     c = (0.0, 0.0, 1.0),
/// ];
/// let atoms = vec![];
/// let crystal = CrystalBuilder::new()
///     .with_lattice(lattice)
///     .with_atoms(&atoms)
///     .build()
///     .unwrap();
/// ```
#[derive(Debug)]
pub struct CrystalBuilder<LatticeState, AtomsState> {
    crystal: Crystal,
    _lattice: std::marker::PhantomData<LatticeState>,
    _atoms: std::marker::PhantomData<AtomsState>,
}

impl CrystalBuilder<LatticeNotSet, AtomsNotSet> {
    #[must_use]
    pub fn new() -> Self {
        CrystalBuilder {
            crystal: Crystal {
                lattice: Lattice::new(
                    [Angstrom(1.0), Angstrom(0.0), Angstrom(0.0)],
                    [Angstrom(0.0), Angstrom(1.0), Angstrom(0.0)],
                    [Angstrom(0.0), Angstrom(0.0), Angstrom(1.0)],
                ),
                positions: vec![],
                kinds: vec![],
            },
            _lattice: std::marker::PhantomData,
            _atoms: std::marker::PhantomData,
        }
    }
}

impl<A> CrystalBuilder<LatticeNotSet, A> {
    // TODO: should Lattice pass as ref?
    #[must_use]
    pub fn with_lattice(self, lattice: Lattice) -> CrystalBuilder<LatticeSet, A> {
        CrystalBuilder {
            crystal: Crystal {
                lattice,
                ..self.crystal
            },
            _lattice: std::marker::PhantomData,
            _atoms: std::marker::PhantomData,
        }
    }
}

impl<L> CrystalBuilder<L, AtomsNotSet> {
    #[must_use]
    pub fn with_atoms(self, atoms: &[Atom]) -> CrystalBuilder<L, AtomsSet> {
        let (positions, kinds) = atoms
            .iter()
            .map(|atom| (atom.position, atom.kind))
            .collect();

        CrystalBuilder {
            crystal: Crystal {
                positions,
                kinds,
                ..self.crystal
            },
            _lattice: std::marker::PhantomData,
            _atoms: std::marker::PhantomData,
        }
    }
}

impl CrystalBuilder<LatticeSet, AtomsSet> {
    /// build and validate the it is a valid crystal.
    /// At the moment only validate the size(positions) == size(numbers)
    pub fn build(self) -> Result<Crystal, CrystalValidateError> {
        // TODO: call validate
        if !self.crystal.positions.len() == self.crystal.kinds.len() {
            return Err(CrystalValidateError {
                message: "crystal valid failed".to_string(),
            });
        }

        Ok(self.crystal)
    }

    // build without runtime validation this is for proc macro which valid in compile time.
    pub(crate) fn build_uncheck(self) -> Crystal {
        self.crystal
    }
}

#[derive(Debug)]
pub struct Atom {
    position: [f64; 3],
    kind: i32,
}

impl Atom {
    // TODO: kind can be more complex, need a type to hold it.
    pub fn new(position: [f64; 3], kind: i32) -> Self {
        Atom { position, kind }
    }
}

// Crystal is the public API so should be align with the real world convention.
// I did not expose the data structure for crystal directly but the builder.
// Internally fileds data structures are private to keep API stable.
// Now I try to align it with moyo's Cell data structure.
//
//
// TODO: not yet generic, but for 3D only at the moment, generalize when doing 2D and 1D
// in prototype, this struct include as much information as possible. may need to be generic.
// use rust's type system to check the problem of frac coordinates and direct coordinates.
#[derive(Debug)]
pub struct Crystal {
    lattice: Lattice,
    positions: Vec<[f64; 3]>,
    kinds: Vec<i32>,
}

impl Crystal {
    #[must_use]
    pub fn builder() -> CrystalBuilder<LatticeNotSet, AtomsNotSet> {
        CrystalBuilder::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn macro_lattice_angstrom() {
        let _ = lattice_angstrom![
            a = (1.0, 0.0, 0.0),
            b = (0.0, 1.0, 0.0),
            c = (0.0, 0.0, 1.0),
        ];
        let _ = lattice_angstrom![(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0),];
        let _ = lattice_angstrom![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],];
        // trailing comma ','
        let _ = lattice_angstrom![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    }

    #[test]
    fn build_crystal_compile_error() {
        let t = trybuild::TestCases::new();
        t.compile_fail("tests/build_crystal/fail_*.rs");
    }
}
