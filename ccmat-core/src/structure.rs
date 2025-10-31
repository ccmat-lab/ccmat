/*
 * base.rs contains or basic structure (crystal and molecule together) utils.
 * notes:
 * Operation safety is guranteed by the type.
 * Use Angstrom as the major internal and default API unit to be consistent with xx/xx.
 * Internally use ``FracCoord`` to represent the position to make sure fractional
 * coords don’t change if lattice changes shape or scale.
 *
 * - Use 'lattice'
 *
 * Compile time build errors include:
 * - when fractional coordinates x not satisfy 0 <= x < 1.0.
 * - multi-set of lattice and sites.
 * - not set lattice or sites.
 *
 * Following errors or runtime validation.
 * - exact duplicate sites (? this might be suitable as compile time error, but how?)
 * - lattice vectors on the same plane
 *
 */

use crate::{
    math::{TransformationMatrix, Vector3},
    symmetry::niggli_reduce,
};

// TODO: naming convention for vars, check IUCr or cif specification
// Give a table to compare in between different popular tools.

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

/// Lattice
/// inner data structure of the struct are private.
/// TODO: this can derive Copy
#[derive(Debug, Clone)]
pub struct Lattice {
    a: Vector3<Angstrom>,
    b: Vector3<Angstrom>,
    c: Vector3<Angstrom>,
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
fn dot(v: &Vector3<f64>, u: &Vector3<f64>) -> f64 {
    v[0] * u[0] + v[1] * u[1] + v[2] * u[2]
}

/// cross product
fn cross(u: &Vector3<f64>, v: &Vector3<f64>) -> Vector3<f64> {
    Vector3::<f64>([
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ])
}

// TODO: Lattice and LatticeReciprocal can be generalized through Unit.

impl Lattice {
    // TODO: how to use type system to validate the row/column definition?
    #[must_use]
    pub fn new(a: Vector3<Angstrom>, b: Vector3<Angstrom>, c: Vector3<Angstrom>) -> Self {
        Lattice { a, b, c }
    }

    #[must_use]
    pub fn a(&self) -> Vector3<Angstrom> {
        self.a
    }

    #[must_use]
    pub fn b(&self) -> Vector3<Angstrom> {
        self.b
    }

    #[must_use]
    pub fn c(&self) -> Vector3<Angstrom> {
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

    #[must_use]
    pub fn volume(&self) -> Volume {
        let (a, b, c) = (self.a.into(), self.b.into(), self.c.into());

        // a⋅(b×c)
        Volume(dot(&a, &cross(&b, &c)))
    }

    #[must_use]
    pub fn reciprocal(&self) -> LatticeReciprocal {
        let (a, b, c) = (self.a.into(), self.b.into(), self.c.into());
        let volume: f64 = Volume(dot(&a, &cross(&b, &c))).into();
        let a_star = 1.0 / volume * (2.0 * std::f64::consts::PI) * cross(&b, &c);
        let b_star = 1.0 / volume * (2.0 * std::f64::consts::PI) * cross(&c, &a);
        let c_star = 1.0 / volume * (2.0 * std::f64::consts::PI) * cross(&a, &b);

        let a_star = Vector3(a_star.map(InvAngstrom::from));
        let b_star = Vector3(b_star.map(InvAngstrom::from));
        let c_star = Vector3(c_star.map(InvAngstrom::from));

        LatticeReciprocal::new(a_star, b_star, c_star)
    }

    /// Find  niggli reduce lattice.
    ///
    /// It using `moyo` to search the niggli reduced lattice, return the reduced lattice and the
    /// operation matrix.
    ///
    /// # Errors
    ///
    /// Error when the search failed which happened if the lattice found is not pass the niggli lattice validation.
    pub fn niggli_reduce(
        &self,
    ) -> Result<(Self, TransformationMatrix), Box<dyn std::error::Error + Send + Sync>> {
        let (a, b, c) = (self.a.into(), self.b.into(), self.c.into());
        let (basis, matrix) = niggli_reduce([a, b, c])?;
        let latt = Lattice::new(basis[0].into(), basis[1].into(), basis[2].into());
        Ok((latt, matrix))
    }

    /// Lattice is represented in the new basis
    ///
    /// (a', b', c') = (a, b, c) * { m00 m01 m02 }
    ///                            { m10 m11 m12 }
    ///                            { m20 m21 m22 }
    ///
    /// a' = m00 * a + m10 * b + m20 * c
    /// b' = m01 * a + m11 * b + m21 * c
    /// c' = m02 * a + m12 * b + m22 * c
    #[must_use]
    pub fn change_basis_by(&self, m: &TransformationMatrix) -> Self {
        let (a, b, c): (Vector3<f64>, Vector3<f64>, Vector3<f64>) =
            (self.a.into(), self.b.into(), self.c.into());

        // dbg!(a, b, c);
        // dbg!(m[0][0], m[1][0], m[2][0]);
        // dbg!(m[0][0] * a);
        // dbg!(m[1][0] * b);
        // dbg!(m[0][0] * a + m[1][0] * b);

        let ap = m[0][0] * a + m[1][0] * b + m[2][0] * c;
        let bp = m[0][1] * a + m[1][1] * b + m[2][1] * c;
        let cp = m[0][2] * a + m[1][2] * b + m[2][2] * c;

        dbg!(ap);
        Self::new(ap.into(), bp.into(), cp.into())
    }
}

// TODO: add lattice_bohr!()
// TODO: impl Display to Lattice for pretty print.

/// Create a [`Lattice`] using vectors expressed in **Ångström** units.
///
/// This macro constructs a [`Lattice`] from three lattice vectors (`a`, `b`, and `c`)
/// expressed as tuples or arrays of three floating-point numbers.
///
/// - Each component is converted to ``Vector3<Angstrom>`` automatically.
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
                $crate::math::Vector3::<$crate::Angstrom>([
                    $crate::Angstrom($x),
                    $crate::Angstrom($y),
                    $crate::Angstrom($z),
                ])
            };
            (($x:expr, $y:expr, $z:expr)) => {
                $crate::math::Vector3::<$crate::Angstrom>([
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
                $crate::math::Vector3::<$crate::Angstrom>([
                    $crate::Angstrom($x),
                    $crate::Angstrom($y),
                    $crate::Angstrom($z),
                ])
            };
            (($x:expr, $y:expr, $z:expr)) => {
                $crate::math::Vector3::<$crate::Angstrom>([
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
    a: Vector3<InvAngstrom>,
    b: Vector3<InvAngstrom>,
    c: Vector3<InvAngstrom>,
}

impl LatticeReciprocal {
    #[must_use]
    pub fn new(
        a_star: Vector3<InvAngstrom>,
        b_star: Vector3<InvAngstrom>,
        c_star: Vector3<InvAngstrom>,
    ) -> Self {
        LatticeReciprocal {
            a: a_star,
            b: b_star,
            c: c_star,
        }
    }

    #[must_use]
    pub fn reciprocal(&self) -> Lattice {
        let (a, b, c) = (self.a.into(), self.b.into(), self.c.into());
        let volume: f64 = Volume(dot(&a, &cross(&b, &c))).into();
        let a_star = 1.0 / volume * (2.0 * std::f64::consts::PI) * cross(&b, &c);
        let b_star = 1.0 / volume * (2.0 * std::f64::consts::PI) * cross(&c, &a);
        let c_star = 1.0 / volume * (2.0 * std::f64::consts::PI) * cross(&a, &b);

        let a_star = Vector3(a_star.map(Angstrom::from));
        let b_star = Vector3(b_star.map(Angstrom::from));
        let c_star = Vector3(c_star.map(Angstrom::from));

        Lattice::new(a_star, b_star, c_star)
    }

    pub fn lattice_params(&self) -> (InvAngstrom, InvAngstrom, InvAngstrom, Rad, Rad, Rad) {
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
            f64::acos(cos_alpha).into(),
            f64::acos(cos_beta).into(),
            f64::acos(cos_gamma).into(),
        )
    }

    pub fn niggli_reduce(
        &self,
    ) -> Result<(Self, TransformationMatrix), Box<dyn std::error::Error + Send + Sync>> {
        let (a, b, c) = (self.a.into(), self.b.into(), self.c.into());
        let (basis, matrix) = niggli_reduce([a, b, c])?;
        let latt = LatticeReciprocal::new(basis[0].into(), basis[1].into(), basis[2].into());
        Ok((latt, matrix))
    }

    /// Lattice is represented in the new basis
    ///
    /// (a', b', c') = (a, b, c) * { m00 m01 m02 }
    ///                            { m10 m11 m12 }
    ///                            { m20 m21 m22 }
    ///
    /// a' = m00 * a + m10 * b + m20 * c
    /// b' = m01 * a + m11 * b + m21 * c
    /// c' = m02 * a + m12 * b + m22 * c
    #[must_use]
    pub fn change_basis_by(&self, m: &TransformationMatrix) -> Self {
        let (a, b, c): (Vector3<f64>, Vector3<f64>, Vector3<f64>) =
            (self.a.into(), self.b.into(), self.c.into());
        let ap = m[0][0] * a + m[1][0] * b + m[2][0] * c;
        let bp = m[0][1] * a + m[1][1] * b + m[2][1] * c;
        let cp = m[0][2] * a + m[1][2] * b + m[2][2] * c;
        Self::new(ap.into(), bp.into(), cp.into())
    }

    #[must_use]
    pub fn a_star(&self) -> Vector3<InvAngstrom> {
        self.a
    }

    #[must_use]
    pub fn b_star(&self) -> Vector3<InvAngstrom> {
        self.b
    }

    #[must_use]
    pub fn c_star(&self) -> Vector3<InvAngstrom> {
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

impl<S> CrystalBuilder<LatticeNotSet, S> {
    // TODO: should Lattice pass as ref?
    #[must_use]
    pub fn with_lattice(self, lattice: &Lattice) -> CrystalBuilder<LatticeSet, S> {
        CrystalBuilder {
            crystal: Crystal {
                // TODO: transfer the ownership instead clone?
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
    // It is also used inside crate where the crystal is known to be valid.
    pub(crate) fn build_uncheck(self) -> Crystal {
        debug_assert!(self.crystal.positions.len() == self.crystal.species.len());

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
pub struct Specie {
    atomic_number: u8,
}

impl Specie {
    fn new(atomic_number: u8) -> Self {
        Specie { atomic_number }
    }

    #[must_use]
    pub fn atomic_number(&self) -> u8 {
        self.atomic_number
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
// Now I try to align it similar to moyo's Cell data structure.
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
        // TODO: is Cow possible?
        self.lattice.clone()
    }

    #[must_use]
    pub fn lattice_reciprocal(&self) -> LatticeReciprocal {
        self.lattice().reciprocal()
    }

    #[must_use]
    pub fn positions(&self) -> Vec<Vector3<FracCoord>> {
        // TODO: avoid clone in readonly?
        self.positions
            .iter()
            .map(|p| Vector3::<FracCoord>(*p))
            .collect()
    }

    #[must_use]
    pub fn species(&self) -> &[Specie] {
        &self.species
    }

    #[must_use]
    pub fn volume(&self) -> Volume {
        self.lattice.volume()
    }
}

// pub fn find_primitive_spglib(
//     crystal: &Crystal,
// ) -> Result<(Crystal, PMatrix, InvPMatrix), Box<dyn std::error::Error + Sync + Send>> {
//     todo!()
// }

#[cfg(test)]
mod tests {
    use crate::{atomic_number, matrix_3x3};

    use super::*;

    #[test]
    fn build_crystal_compile_error() {
        let t = trybuild::TestCases::new();
        t.compile_fail("tests/build_crystal/^fail_*.rs");
    }

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
    fn reciprocal() {
        let lattice = lattice_angstrom![
            // no orthogonal cell
            a = (2.0, 0.5, 0.0),
            b = (0.0, 3.0, 0.5),
            c = (0.5, 0.0, 4.0),
        ];

        let latt2 = lattice.reciprocal().reciprocal();

        // TODO: an assert with assert_approx! util
        dbg!(latt2, lattice);
    }

    #[test]
    fn crystal_volume() {
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
    }

    #[test]
    fn change_basis_by() {
        let lattice = lattice_angstrom![
            // no orthogonal cell
            a = (2.0, 0.5, 0.0),
            b = (0.0, 3.0, 0.5),
            c = (0.5, 0.0, 4.0),
        ];

        let tmatrix = matrix_3x3![
            1 0 0;
            0 1 0;
            0 0 1;
        ];

        dbg!(lattice.change_basis_by(&tmatrix));

        // let tmatrix = matrix_3x3![
        //     sqrt(3)/2  -1/2     0;
        //     1/2   sqrt(3)/2     0;
        //     0             0     1;
        // ];
    }
}
