mod atomic;

use std::{
    borrow::Cow,
    ops::{Deref, Mul},
};

use moyo::{
    self,
    data::{arithmetic_crystal_class_entry, hall_symbol_entry},
    MoyoDataset,
};
use nalgebra::Vector3;

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

impl From<Angstrom> for f64 {
    fn from(value: Angstrom) -> Self {
        value.0
    }
}

impl From<f64> for Angstrom {
    fn from(value: f64) -> Self {
        Angstrom(value)
    }
}

/// Unit the inverse Angstrom 1/Å for vectors in reciprocal space.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InvAngstrom(pub f64);

impl From<InvAngstrom> for f64 {
    fn from(value: InvAngstrom) -> Self {
        value.0
    }
}

impl From<f64> for InvAngstrom {
    fn from(value: f64) -> Self {
        InvAngstrom(value)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bohr(pub f64);

impl From<Bohr> for f64 {
    fn from(value: Bohr) -> Self {
        value.0
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FracCoord(pub f64);

impl std::fmt::Display for FracCoord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:15.9}", self.0)
    }
}

impl From<FracCoord> for f64 {
    fn from(value: FracCoord) -> Self {
        value.0
    }
}

impl From<f64> for FracCoord {
    fn from(value: f64) -> Self {
        FracCoord(value)
    }
}

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

// #[macro_export]
// macro_rules! atomic_number {
//     (I)  => { 53 };
//     (Cu) => { 29 };
// }

/// macro to set the sites
///
/// # Examples
///
/// ```
/// use ccmat_core::sites_frac_coord;
///
/// let _ = sites_frac_coord![
///     (0.0, 0.0, 0.0), 8;
///     (0.0, 0.0, 0.5), 8;
/// ];
/// ```
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
/// TODO: this can derive Copy
#[derive(Debug, Clone)]
pub struct Lattice {
    a: Vec3<Angstrom>,
    b: Vec3<Angstrom>,
    c: Vec3<Angstrom>,
}

/// f64 wrapper for radians
#[derive(Debug, Copy, Clone)]
pub struct Rad(f64);

impl From<Rad> for f64 {
    fn from(value: Rad) -> Self {
        value.0
    }
}

impl From<f64> for Rad {
    fn from(value: f64) -> Self {
        Rad(value)
    }
}

// TODO: I should have a proc-macro for impl all such from f64 traits

/// f64 wrapper for value with unit of volume (Angstrom ^ 2)
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Volume(f64);

impl From<Volume> for f64 {
    fn from(value: Volume) -> Self {
        value.0
    }
}

impl From<f64> for Volume {
    fn from(value: f64) -> Self {
        Volume(value)
    }
}

/// dot product
fn dot(v: &Vec3<f64>, u: &Vec3<f64>) -> f64 {
    v[0] * u[0] + v[1] * u[1] + v[2] + u[2]
}

/// cross product
fn cross(u: &Vec3<f64>, v: &Vec3<f64>) -> Vec3<f64> {
    Vec3::<f64>([
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ])
}

#[derive(Debug, Copy, Clone)]
pub struct Vec3<T>(pub [T; 3]);

impl<T> Deref for Vec3<T> {
    type Target = [T; 3];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T> std::ops::Index<usize> for Vec3<T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl Mul<Vec3<f64>> for f64 {
    type Output = Vec3<f64>;

    fn mul(self, rhs: Vec3<f64>) -> Self::Output {
        let v = rhs.0.map(f64::from);
        Vec3::<f64>(v)
    }
}

impl From<Vec3<Angstrom>> for Vec3<f64> {
    fn from(v: Vec3<Angstrom>) -> Self {
        Vec3::<f64>(v.map(f64::from))
    }
}

impl Lattice {
    // TODO: how to use type system to validate the row/column definition?
    #[must_use]
    pub fn new(a: Vec3<Angstrom>, b: Vec3<Angstrom>, c: Vec3<Angstrom>) -> Self {
        Lattice { a, b, c }
    }

    #[must_use]
    pub fn a(&self) -> Vec3<Angstrom> {
        self.a
    }

    #[must_use]
    pub fn b(&self) -> Vec3<Angstrom> {
        self.b
    }

    #[must_use]
    pub fn c(&self) -> Vec3<Angstrom> {
        self.c
    }

    pub fn lattice_params(&self) -> (Angstrom, Angstrom, Angstrom, Rad, Rad, Rad) {
        let va = self.a.map(f64::from);
        let length_a = f64::sqrt(va[0] * va[0] + va[1] * va[1] + va[2] * va[2]);

        let vb = self.b.map(f64::from);
        let length_b = f64::sqrt(vb[0] * vb[0] + vb[1] * vb[1] + vb[2] * vb[2]);

        let vc = self.c.map(f64::from);
        let length_c = f64::sqrt(vc[0] * vc[0] + vc[1] * vc[1] + vc[2] * vc[2]);

        let cos_alpha = (vb[0] * vc[0] + vb[1] * vc[1] + vb[2] * vc[2]) / (length_b * length_c);
        let cos_beta = (va[0] * vc[0] + va[1] * vc[1] + va[2] * vc[2]) / (length_a * length_c);
        let cos_gamma = (va[0] * vb[0] + va[1] * vb[1] + va[2] * vb[2]) / (length_a * length_b);

        (
            length_a.into(),
            length_b.into(),
            length_c.into(),
            cos_alpha.acos().into(),
            cos_beta.acos().into(),
            cos_gamma.acos().into(),
        )
    }

    #[must_use]
    pub fn length_a(&self) -> Angstrom {
        self.lattice_params().0
    }

    #[must_use]
    pub fn length_b(&self) -> Angstrom {
        self.lattice_params().1
    }

    #[must_use]
    pub fn length_c(&self) -> Angstrom {
        self.lattice_params().2
    }

    #[must_use]
    pub fn rad_alpha(&self) -> Rad {
        self.lattice_params().3
    }

    #[must_use]
    pub fn rad_beta(&self) -> Rad {
        self.lattice_params().4
    }

    #[must_use]
    pub fn rad_gamma(&self) -> Rad {
        self.lattice_params().5
    }

    pub fn volume(&self) -> Volume {
        let (a, b, c) = (self.a.into(), self.b.into(), self.c.into());

        // a⋅(b×c)
        Volume(dot(&a, &cross(&b, &c)))
    }

    #[must_use]
    pub fn reciprocal_lattice(&self) -> LatticeReciprocal {
        let (a, b, c) = (self.a.into(), self.b.into(), self.c.into());
        let volume: f64 = Volume(dot(&a, &cross(&b, &c))).into();
        let a_star = 1.0 / volume * (2.0 * std::f64::consts::PI) * cross(&b, &c);
        let b_star = 1.0 / volume * (2.0 * std::f64::consts::PI) * cross(&c, &a);
        let c_star = 1.0 / volume * (2.0 * std::f64::consts::PI) * cross(&a, &b);

        let a_star = a_star.map(InvAngstrom::from);
        let b_star = b_star.map(InvAngstrom::from);
        let c_star = c_star.map(InvAngstrom::from);

        LatticeReciprocal::new(a_star, b_star, c_star)
    }
}

// TODO: add lattice_bohr!()
// TODO: impl Display to Lattice for pretty print.

/// Create a [`Lattice`] using vectors expressed in **Ångström** units.
///
/// This macro constructs a [`Lattice`] from three lattice vectors (`a`, `b`, and `c`)
/// expressed as tuples or arrays of three floating-point numbers.
///
/// - Each component is converted to ``Vec3<Angstrom>`` automatically.
/// - Both `(x, y, z)` and `[x, y, z]` tuple/array syntax are supported for vector.
/// - Trailing commas are optional.
///
/// It supports both **named** and **positional** forms:
///
/// - **Named form** (explicit `a=`, `b=`, `c=`):
///   ```
///   use ccmat_core::lattice_angstrom;
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
///   use ccmat_core::{lattice_angstrom, Lattice, Angstrom};
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
///
/// It is also okay to use '[]' instead of '()' for each vector.
///
/// ```
/// use ccmat_core::{lattice_angstrom, Lattice, Angstrom};
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
        macro_rules! __vec3_angstrom {
            ([$x:expr, $y:expr, $z:expr]) => {
                $crate::Vec3::<$crate::Angstrom>([
                    $crate::Angstrom($x),
                    $crate::Angstrom($y),
                    $crate::Angstrom($z),
                ])
            };
            (($x:expr, $y:expr, $z:expr)) => {
                $crate::Vec3::<$crate::Angstrom>([
                    $crate::Angstrom($x),
                    $crate::Angstrom($y),
                    $crate::Angstrom($z),
                ])
            };
        }

        let lattice = $crate::Lattice::new(
            __vec3_angstrom!($a),
            __vec3_angstrom!($b),
            __vec3_angstrom!($c),
        );
        lattice
    }};
    (
        $a:tt,
        $b:tt,
        $c:tt $(,)?
    ) => {{
        macro_rules! __vec3_angstrom {
            ([$x:expr, $y:expr, $z:expr]) => {
                $crate::Vec3::<$crate::Angstrom>([
                    $crate::Angstrom($x),
                    $crate::Angstrom($y),
                    $crate::Angstrom($z),
                ])
            };
            (($x:expr, $y:expr, $z:expr)) => {
                $crate::Vec3::<$crate::Angstrom>([
                    $crate::Angstrom($x),
                    $crate::Angstrom($y),
                    $crate::Angstrom($z),
                ])
            };
        }
        let lattice = $crate::Lattice::new(
            __vec3_angstrom!($a),
            __vec3_angstrom!($b),
            __vec3_angstrom!($c),
        );
        lattice
    }};
}

pub struct LatticeReciprocal {
    // internal use a not a_star, but the API is a_star to make it very explicit.
    a: [InvAngstrom; 3],
    b: [InvAngstrom; 3],
    c: [InvAngstrom; 3],
}

impl LatticeReciprocal {
    #[must_use]
    pub fn new(
        a_star: [InvAngstrom; 3],
        b_star: [InvAngstrom; 3],
        c_star: [InvAngstrom; 3],
    ) -> Self {
        LatticeReciprocal {
            a: a_star,
            b: b_star,
            c: c_star,
        }
    }

    #[must_use]
    pub fn a_star(&self) -> [InvAngstrom; 3] {
        self.a
    }

    #[must_use]
    pub fn b_star(&self) -> [InvAngstrom; 3] {
        self.b
    }

    #[must_use]
    pub fn c_star(&self) -> [InvAngstrom; 3] {
        self.c
    }
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
/// use ccmat_core::*;
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

    #[must_use]
    pub fn lattice(&self) -> Lattice {
        // TODO: Cow??
        self.lattice.clone()
    }

    #[must_use]
    pub fn volume(&self) -> Volume {
        self.lattice.volume()
    }
}

// The reference version not provide, since the moyo::Cell only used internally thus
// assume no need to hold its ownership.
impl TryFrom<moyo::base::Cell> for Crystal {
    type Error = Box<dyn std::error::Error + Sync + Send>;

    fn try_from(s: moyo::base::Cell) -> Result<Self, Self::Error> {
        // TODO: suggest API to get every basis in ccmat-moyo
        let a = s.lattice.basis.column(0);
        let b = s.lattice.basis.column(1);
        let c = s.lattice.basis.column(2);

        let lattice = lattice_angstrom![
            a = (a[0], a[1], a[2]),
            b = (b[0], b[1], b[2]),
            c = (c[0], c[1], c[2]),
        ];

        let positions = s.positions;
        let numbers = s.numbers;

        // Length are guranteed to be the same because the moyo::Cell is constructed from ccmat
        // use `into_iter` to move and avoid allocation.
        let sites: Vec<Site> = positions
            .into_iter()
            .zip(numbers)
            .map(|(pos, num)| {
                let pos: [FracCoord; 3] = [pos[0].into(), pos[1].into(), pos[2].into()];
                let num: u8 = num.try_into().expect("atomic number not in 0..128");
                Site::new(pos, num)
            })
            .collect();

        let c = CrystalBuilder::new()
            .with_lattice(&lattice)
            .with_sites(&sites)
            .build()?;
        Ok(c)
    }
}

// TODO: pr to moyo ask for a constructor API, experiment it in ccmat-moyo first.
impl From<&Crystal> for moyo::base::Cell {
    fn from(s: &Crystal) -> Self {
        let a = s.lattice.a().map(f64::from);
        let b = s.lattice.b().map(f64::from);
        let c = s.lattice.c().map(f64::from);
        let lattice = moyo::base::Lattice::from_basis([a, b, c]);

        // TODO: moyo need an api or macro to create positions.
        let positions = s
            .positions
            .iter()
            // TODO: ccmat-moyo, add an API `position!` for Vector3.
            .map(|p| Vector3::new(p[0].into(), p[1].into(), p[2].into()))
            .collect();

        let numbers = s.species.iter().map(|s| s.atomic_number.into()).collect();

        moyo::base::Cell::new(lattice, positions, numbers)
    }
}

// If willing to give the ownership.
impl From<Crystal> for moyo::base::Cell {
    fn from(s: Crystal) -> Self {
        (&s).into()
    }
}

#[allow(non_camel_case_types)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BravaisClass {
    // Triclinic
    aP,
    // Monoclinic
    mP,
    mC,
    // Orthorhombic
    oP,
    oS,
    oF,
    oI,
    // Tetragonal
    tP,
    tI,
    // Rhombohedral
    hR,
    // Hexagonal
    hP,
    // Cubic
    cP,
    cF,
    cI,
}

impl From<moyo::data::BravaisClass> for BravaisClass {
    fn from(bv: moyo::data::BravaisClass) -> Self {
        match bv {
            moyo::data::BravaisClass::aP => BravaisClass::aP,
            moyo::data::BravaisClass::mP => BravaisClass::mP,
            moyo::data::BravaisClass::mC => BravaisClass::mC,
            moyo::data::BravaisClass::oP => BravaisClass::oP,
            moyo::data::BravaisClass::oS => BravaisClass::oS,
            moyo::data::BravaisClass::oF => BravaisClass::oF,
            moyo::data::BravaisClass::oI => BravaisClass::oI,
            moyo::data::BravaisClass::tP => BravaisClass::tP,
            moyo::data::BravaisClass::tI => BravaisClass::tI,
            moyo::data::BravaisClass::hR => BravaisClass::hR,
            moyo::data::BravaisClass::hP => BravaisClass::hP,
            moyo::data::BravaisClass::cP => BravaisClass::cP,
            moyo::data::BravaisClass::cF => BravaisClass::cF,
            moyo::data::BravaisClass::cI => BravaisClass::cI,
        }
    }
}

impl std::fmt::Display for BravaisClass {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            // Triclinic
            BravaisClass::aP => "aP",
            // Monoclinic
            BravaisClass::mP => "mP",
            BravaisClass::mC => "mC",
            // Orthorhombic
            BravaisClass::oP => "oP",
            BravaisClass::oS => "oS",
            BravaisClass::oF => "oF",
            BravaisClass::oI => "oI",
            // Tetragonal
            BravaisClass::tP => "tP",
            BravaisClass::tI => "tI",
            // Rhombohedral
            BravaisClass::hR => "hR",
            // Hexagonal
            BravaisClass::hP => "hP",
            // Cubic
            BravaisClass::cP => "cP",
            BravaisClass::cF => "cF",
            BravaisClass::cI => "cI",
        };
        write!(f, "{s}")
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum Centering {
    P, // Primitive
    A, // A-face centered
    B, // B-face centered
    C, // C-face centered
    I, // Body centered
    R, // Rhombohedral (obverse setting)
    F, // Face centered
}

impl From<moyo::data::Centering> for Centering {
    fn from(c: moyo::data::Centering) -> Self {
        match c {
            moyo::data::Centering::P => Centering::P,
            moyo::data::Centering::A => Centering::A,
            moyo::data::Centering::B => Centering::B,
            moyo::data::Centering::C => Centering::C,
            moyo::data::Centering::I => Centering::I,
            moyo::data::Centering::R => Centering::R,
            moyo::data::Centering::F => Centering::F,
        }
    }
}

/// Wrapper of `MoyoDataset` with handy APIs.
pub struct SymmetryInfo {
    inner: MoyoDataset,
}

impl SymmetryInfo {
    /// Space group number (1-230)
    ///
    /// # Panics
    /// Panic if space group number return from moyo is negative, i32 -> u32 fail, should be a bug
    /// then.
    #[must_use]
    pub fn spg_number(&self) -> u32 {
        self.inner
            .number
            .try_into()
            .expect("spage group number not in 1..=230")
    }

    /// Hall symbol number (1-530).
    ///
    /// # Panics
    /// Panic if hall number return from moyo is negative, i32 -> u32 fail, should be a bug
    /// then.
    #[must_use]
    pub fn hall_number(&self) -> u32 {
        self.inner
            .hall_number
            .try_into()
            .expect("spage group number not in 1..=230")
    }

    /// Bravais class
    ///
    /// # Panics
    /// When moyo failed to get the hall symbol from the `hall_number`, shouldn't happened in ccmat
    /// since the `hall_number` is computed by moyo, if panic it is a bug.
    #[must_use]
    pub fn bravais_class(&self) -> BravaisClass {
        let hall_number = self.inner.hall_number;
        let hall_symbol =
            hall_symbol_entry(hall_number).expect("unable to get hall symbol from hall_number");
        let arithmetic_entry =
            arithmetic_crystal_class_entry(hall_symbol.arithmetic_number).unwrap();

        arithmetic_entry.bravais_class.into()
    }

    /// Centering
    ///
    /// # Panics
    /// When moyo failed to get the hall symbol from the `hall_number`, shouldn't happened in ccmat
    /// since the `hall_number` is computed by moyo, if panic it is a bug.
    #[must_use]
    pub fn centring(&self) -> Centering {
        let hall_number = self.inner.hall_number;
        let hall_symbol =
            hall_symbol_entry(hall_number).expect("unable to get hall symbol from hall_number");
        hall_symbol.centering.into()
    }

    /// Hall symbol
    ///
    /// # Panics
    /// When moyo failed to get the hall symbol from the `hall_number`, shouldn't happened in ccmat
    /// since the `hall_number` is computed by moyo, if panic it is a bug.
    #[must_use]
    pub fn hall_symbol(&self) -> Cow<'_, str> {
        let hall_number = self.inner.hall_number;
        let hall_symbol =
            hall_symbol_entry(hall_number).expect("unable to get hall symbol from hall_number");
        Cow::Borrowed(hall_symbol.hall_symbol)
    }

    /// Check if contain inversion symmetry.
    #[must_use]
    pub fn has_inversion(&self) -> bool {
        let hall_symbol = self.hall_symbol();
        hall_symbol.starts_with('-')
    }

    /// Stdandard Cell
    ///
    /// # Errors
    /// ??
    pub fn std_cell(&self) -> Result<Crystal, Box<dyn std::error::Error + Send + Sync>> {
        self.inner.std_cell.clone().try_into()
    }
}

// pub fn find_primitive_spglib(
//     crystal: &Crystal,
// ) -> Result<(Crystal, PMatrix, InvPMatrix), Box<dyn std::error::Error + Sync + Send>> {
//     todo!()
// }

// TODO: move to ccmat-moyo and add wrapper for MagneticMoyoDataset, and PR to moyo to integrate.
// It can be generic over Mag/nonmag dataset.

/// analyze symmetry
///
/// # Errors
/// ???
// TODO: thiserror
pub fn analyze_symmetry(
    crystal: &Crystal,
    symprec: f64,
) -> Result<SymmetryInfo, Box<dyn std::error::Error + Send + Sync>> {
    let cell: moyo::base::Cell = crystal.into();
    let sym_info = MoyoDataset::with_default(&cell, symprec)?;
    let sym_info = SymmetryInfo { inner: sym_info };
    Ok(sym_info)
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
    fn moyo_spg() {
        const X_4F: f64 = 0.3046;

        let lattice = lattice_angstrom![
            a = (4.603, 0.0, 0.0),
            b = (0.0, 4.603, 0.0),
            c = (0.0, 0.0, 4.603),
        ];
        let sites = sites_frac_coord![
            (0.0, 0.0, 0.0), atomic_number!(Ti);               // Ti(2a)
            (0.5, 0.5, 0.5), atomic_number!(Ti);               // Ti(2a)
            (X_4F, X_4F, 0.0), atomic_number!(O);              // O(4f)
            (1.0 - X_4F, 1.0 - X_4F, 0.0), atomic_number!(O);  // O(4f)
            (-X_4F + 0.5, X_4F + 0.5, 0.5), atomic_number!(O); // O(4f)
            (X_4F + 0.5, -X_4F + 0.5, 0.5), atomic_number!(O); // O(4f)
        ];
        let crystal = CrystalBuilder::new()
            .with_lattice(&lattice)
            .with_sites(&sites)
            .build()
            .unwrap();

        assert_eq!(crystal.volume(), Volume(4.603 * 4.603 * 4.603));

        let syminfo = analyze_symmetry(&crystal, 1e-4).unwrap();
        assert_eq!(syminfo.spg_number(), 136);
    }

    #[test]
    fn build_crystal_compile_error() {
        let t = trybuild::TestCases::new();
        t.compile_fail("tests/build_crystal/^fail_*.rs");
    }
}
