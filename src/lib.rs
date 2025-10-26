/// notes:
/// Operation safety is guranteed by the type.
/// Use Angstrom as the major internal and default API unit to be consistent with xx/xx.
/// Internally use ``FracCoord`` to represent the position to make sure fractional
/// coords don’t change if lattice changes shape or scale.
///
/// - Use 'lattice'
///
/// Compile time build errors include:
/// - when fractional coordinates x not satisfy 0 <= x < 1.0.
/// - multi-set of lattice and sites.
/// - not set lattice or sites.
///
/// Following errors or runtime validation.
/// - exact duplicate sites (? this might be suitable as compile time error, but how?)
/// - lattice vectors on the same plane

// TODO: naming convention for vars, check IUCr or cif specification
// Give a table to compare in between different popular tools.
//
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Angstrom(pub f64);

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bohr(pub f64);

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FracCoord(pub f64);

/// `frac!` macro to create `FracCoord` and validate the value is in between [0.0, 1.0)
/// at compile time.
#[macro_export]
macro_rules! frac {
    ($x:expr) => {{
        let frac_coord = $crate::FracCoord($x);
        const {
            assert!(
                (0.0 <= $x && $x < 1.0),
                "invalid fractional coordinate: must satisfy 0.0 <= x < 1.0"
            );
        }
        frac_coord
    }};
}

#[macro_export]
macro_rules! sites_frac_coord {
    () => {
        Vec::new()
    };
    ( $(
        ($x:expr,$y:expr,$z:expr), $kind:expr
      );+ $(;)?
    ) => {{
        let sites = vec![
            $(
                $crate::Site::new(
                    [
                        $crate::frac!($x),
                        $crate::frac!($y),
                        $crate::frac!($z),
                    ],
                    $kind,
                )
            ),+
        ];
        sites
    }};
}

/// Lattice
/// inner data structure of the struct are private.
#[derive(Debug, Clone)]
pub struct Lattice {
    a: [Angstrom; 3],
    b: [Angstrom; 3],
    c: [Angstrom; 3],
}

impl Lattice {
    // TODO: how to use type system to validate the row/column definition or unit?
    #[must_use]
    pub fn new(a: [Angstrom; 3], b: [Angstrom; 3], c: [Angstrom; 3]) -> Self {
        Lattice { a, b, c }
    }

    #[must_use]
    pub fn a(&self) -> [Angstrom; 3] {
        self.a
    }

    #[must_use]
    pub fn b(&self) -> [Angstrom; 3] {
        self.b
    }

    #[must_use]
    pub fn c(&self) -> [Angstrom; 3] {
        self.c
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
///   use ccmat::lattice_angstrom;
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
///   use ccmat::{lattice_angstrom, Lattice, Angstrom};
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
/// use ccmat::{lattice_angstrom, Lattice, Angstrom};
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
pub struct SitesSet;
pub struct SitesNotSet;

/// use builder pattern so the validation is runtime
///
/// # Example
/// ```
/// use ccmat::*;
///
/// let lattice = lattice_angstrom![
///     a = (1.0, 0.0, 0.0),
///     b = (0.0, 1.0, 0.0),
///     c = (0.0, 0.0, 1.0),
/// ];
/// let sites = vec![];
/// let crystal = CrystalBuilder::new()
///     .with_lattice(&lattice)
///     .with_sites(&sites)
///     .build()
///     .unwrap();
/// ```
#[derive(Debug)]
pub struct CrystalBuilder<LatticeSetState, SiteSetState> {
    crystal: Crystal,
    _lattice: std::marker::PhantomData<LatticeSetState>,
    _sites: std::marker::PhantomData<SiteSetState>,
}

impl Default for CrystalBuilder<LatticeNotSet, SitesNotSet> {
    fn default() -> Self {
        Self {
            crystal: Crystal {
                lattice: lattice_angstrom!([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0],),
                positions: vec![],
                species: vec![],
            },
            _lattice: std::marker::PhantomData,
            _sites: std::marker::PhantomData,
        }
    }
}

impl CrystalBuilder<LatticeNotSet, SitesNotSet> {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }
}

impl<A> CrystalBuilder<LatticeNotSet, A> {
    // TODO: should Lattice pass as ref?
    #[must_use]
    pub fn with_lattice(self, lattice: &Lattice) -> CrystalBuilder<LatticeSet, A> {
        CrystalBuilder {
            crystal: Crystal {
                lattice: lattice.clone(),
                ..self.crystal
            },
            _lattice: std::marker::PhantomData,
            _sites: std::marker::PhantomData,
        }
    }
}

impl<L> CrystalBuilder<L, SitesNotSet> {
    #[must_use]
    pub fn with_sites(self, sites: &[Site]) -> CrystalBuilder<L, SitesSet> {
        let (positions, species) = sites
            .iter()
            .map(|atom| (atom.position, atom.specie.clone()))
            .collect();

        CrystalBuilder {
            crystal: Crystal {
                positions,
                species,
                ..self.crystal
            },
            _lattice: std::marker::PhantomData,
            _sites: std::marker::PhantomData,
        }
    }
}

impl CrystalBuilder<LatticeSet, SitesSet> {
    fn validate(&self) -> Result<(), CrystalValidateError> {
        // TODO: call validate
        if !self.crystal.positions.len() == self.crystal.species.len() {
            return Err(CrystalValidateError {
                message: "crystal valid failed".to_string(),
            });
        }
        Ok(())
    }

    // build without runtime validation this is for proc macro which valid in compile time.
    // uncheck logic errors
    fn build_uncheck(self) -> Crystal {
        self.crystal
    }

    /// build and validate the it is a valid crystal.
    /// At the moment only validate the size(positions) == size(numbers)
    ///
    /// # Errors
    /// ??
    pub fn build(self) -> Result<Crystal, CrystalValidateError> {
        self.validate()?;

        let crystal = self.build_uncheck();

        Ok(crystal)
    }
}

// TODO: partial occupation on sites
#[derive(Debug, Clone)]
struct Specie {
    atomic_number: u8,
}

impl Specie {
    fn new(atomic_number: u8) -> Self {
        Specie { atomic_number }
    }
}

#[derive(Debug)]
pub struct Site {
    position: [FracCoord; 3],
    specie: Specie,
}

impl Site {
    #[must_use]
    pub fn new(position: [FracCoord; 3], atomic_number: u8) -> Self {
        Site {
            position,
            specie: Specie::new(atomic_number),
        }
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
    positions: Vec<[FracCoord; 3]>,
    species: Vec<Specie>,
}

impl Crystal {
    #[must_use]
    pub fn builder() -> CrystalBuilder<LatticeNotSet, SitesNotSet> {
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
    fn macro_sites_frac() {
        let _: Vec<Site> = sites_frac_coord![];
        let _ = sites_frac_coord![
            (0.0, 0.0, 0.0), 8;
            (0.0, 0.0, 0.5), 8;
        ];
    }

    #[test]
    fn build_crystal_compile_error() {
        let t = trybuild::TestCases::new();
        t.compile_fail("tests/build_crystal/^fail_*.rs");
    }
}
