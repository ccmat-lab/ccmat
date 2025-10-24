// TODO: naming convention for vars, check IUCr or cif specification
// Give a table to compare in between different popular tools.

#[derive(Debug, Default)]
pub enum Unit {
    #[default]
    Metric,
    Atomic,
    Relative, // XXX: naming?
}

#[derive(Debug, Default)]
pub enum CoordinateSystem {
    #[default]
    Cardition,
    Relative, // XXX: naming?
}

#[derive(Debug, Default)]
pub struct Lattice {
    basis: [[f64; 3]; 3],
    unit: Unit,
    coord_system: CoordinateSystem,
}

impl Lattice {
    // TODO: how to use type system to validate the row/column definition or unit?
    pub fn new(basis: [[f64; 3]; 3], unit: Unit, coord_system: CoordinateSystem) -> Self {
        Lattice {
            basis,
            unit,
            coord_system,
        }
    }
}

#[derive(Debug)]
struct CrystalValidateError {
    message: String,
}

impl std::fmt::Display for CrystalValidateError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl std::error::Error for CrystalValidateError {}

// use builder pattern so the validation is runtime
// To make it a compile time check, I can use the proc macro (build/build_unchcek API).
#[derive(Debug)]
pub struct CrystalBuilder {
    crystal: Crystal,
}

impl CrystalBuilder {
    #[must_use]
    pub fn new() -> Self {
        CrystalBuilder {
            crystal: Crystal {
                lattice: Lattice::default(),
                positions: vec![],
                kinds: vec![],
            },
        }
    }

    // TODO: should Lattice pass as ref?
    #[must_use]
    pub fn with_lattice(self, lattice: Lattice) -> Self {
        CrystalBuilder {
            crystal: Crystal {
                lattice,
                ..self.crystal
            },
        }
    }

    #[must_use]
    pub fn with_atoms(self, atoms: &[Atom]) -> Self {
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
        }
    }

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
    fn new(position: [f64; 3], kind: i32) -> Self {
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
#[derive(Debug, Default)]
pub struct Crystal {
    lattice: Lattice,
    positions: Vec<[f64; 3]>,
    kinds: Vec<i32>,
}

impl Crystal {
    #[must_use]
    pub fn builder() -> CrystalBuilder {
        CrystalBuilder::new()
    }
}

fn main() {
    let lattice = Lattice::default();
    let atoms = vec![];
    let crystal = CrystalBuilder::new()
        .with_lattice(lattice)
        .with_atoms(&atoms)
        .build();

    dbg!(crystal);

    // this should be a compiler error since it in not fully initialized
    // use state to annotate
    let lattice = Lattice::default();
    let crystal = CrystalBuilder::new().with_lattice(lattice).build();

    dbg!(crystal);
}
