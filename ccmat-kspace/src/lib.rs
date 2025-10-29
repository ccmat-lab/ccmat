mod path;

use log::warn;
use ccmat_core::{analyze_symmetry, BravaisClass, Crystal};

use crate::path::KpathInfo;

fn find_primitive_hpkot(
    cell_std: &Crystal,
    symprec: f64,
) -> Result<moyo::base::Cell, Box<dyn std::error::Error + Send + Sync>> {
    todo!()
}

#[allow(non_camel_case_types)]
#[derive(Debug)]
pub enum ExtBravaisClass {
    // Triclinic
    aP1, // reserved for aP2 + aP3, ref: hpkot paper (Table 94).
    aP2,
    aP3,
    // Monoclinic
    mP1,
    mC1,
    mC2,
    mC3,
    // Orthorhombic
    oP1,
    oA1,
    oA2,
    oC1,
    oC2,
    oF1,
    oF2,
    oF3,
    oI1,
    oI2,
    oI3,
    // Tetragonal
    tP1,
    tI1,
    tI2,
    // Rhombohedral
    hR1,
    hR2,
    // Hexagonal
    hP1,
    hP2,
    // Cubic
    cP1,
    cP2,
    cF1,
    cF2,
    cI1,
}

/// # Errors
/// ??
///
/// # Panics
/// ??
#[allow(clippy::too_many_lines)]
pub fn find_path(
    crystal: &Crystal,
    symprec: f64,
    threshold: f64,
) -> Result<(&'static KpathInfo, Crystal), Box<dyn std::error::Error + Send + Sync>> {
    let syminfo = analyze_symmetry(crystal, symprec)?;
    let cell_std = syminfo.std_cell()?;
    let spg_number = syminfo.spg_number();

    let crystal_priv: Crystal = find_primitive_hpkot(&cell_std, symprec)?.try_into()?;
    let lattice_params = crystal_priv.lattice().lattice_params();
    let (a, b, c, alpha, beta, gamma) = lattice_params;
    let a: f64 = a.into();
    let b: f64 = b.into();
    let c: f64 = c.into();
    let alpha: f64 = alpha.into();
    let beta: f64 = beta.into();
    let gamma: f64 = gamma.into();

    let ext_bravais = match syminfo.bravais_class() {
        BravaisClass::aP => todo!(),
        BravaisClass::mP => ExtBravaisClass::mP1,
        BravaisClass::mC => {
            let cosbeta = f64::cos(beta);
            if f64::abs(b - a * f64::sqrt(1.0 - cosbeta * cosbeta)) < threshold {
                warn!("mC lattice, but b ~ a*sin(beta)");
            }

            if b < a * f64::sqrt(1.0 - cosbeta * cosbeta) {
                ExtBravaisClass::mC1
            } else {
                if f64::abs(-a * cosbeta / c + (a * a) * (1.0 - cosbeta * cosbeta) / (b * b) - 1.0)
                    < threshold
                {
                    warn!("mC lattice, but -a*cos(beta)/c + a^2*sin(beta)^2/b^2 ~ 1");
                }

                if -a * cosbeta / c + (a * a) * (1.0 - cosbeta * cosbeta) / (b * b) < 1.0 {
                    ExtBravaisClass::mC2
                } else {
                    ExtBravaisClass::mC3
                }
            }
        }
        BravaisClass::oP => ExtBravaisClass::oP1,
        BravaisClass::oS => match spg_number {
            // oA
            x if (38..=41).contains(&x) => {
                if f64::abs(b - c) < threshold {
                    warn!("oA lattice, but b ~ c");
                }
                if b < c {
                    ExtBravaisClass::oA1
                } else {
                    ExtBravaisClass::oA2
                }
            }
            // oC
            x if (20..=21).contains(&x) || (35..=37).contains(&x) || (63..=68).contains(&x) => {
                if f64::abs(b - a) < threshold {
                    warn!("oC lattice, but a ~ b");
                }
                if a < b {
                    ExtBravaisClass::oC1
                } else {
                    ExtBravaisClass::oC2
                }
            }
            _ => unreachable!("oS bravais lattice spacegroup number in wrong range"),
        },
        BravaisClass::oF => {
            if f64::abs((1.0 / a * a) - ((1.0 / b * b) + (1.0 / c * c))) < threshold {
                warn!("oF lattice, but 1/a^2 ~ 1/b^2 + 1/c^2");
            }
            if f64::abs((1.0 / c * c) - ((1.0 / a * a) + (1.0 / b * b))) < threshold {
                warn!("oF lattice, but 1/c^2 ~ 1/a^2 + 1/b^2");
            }
            if 1.0 / a * a > (1.0 / b * b) + (1.0 / c * c) {
                ExtBravaisClass::oF1
            } else if 1.0 / c * c > (1.0 / a * a) + (1.0 / b * b) {
                ExtBravaisClass::oF2
            } else {
                ExtBravaisClass::oF3
            }
        }
        BravaisClass::oI => {
            #[derive(Debug)]
            enum Face {
                A, // oI2
                B, // oI3
                C, // oI1
            }
            impl std::fmt::Display for Face {
                fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                    match self {
                        Face::A => write!(f, "A"),
                        Face::B => write!(f, "B"),
                        Face::C => write!(f, "C"),
                    }
                }
            }

            let mut vec = vec![(a, Face::A), (b, Face::B), (c, Face::C)];
            vec.sort_by(|x, y| {
                y.0.partial_cmp(&x.0)
                    .expect("lattice length compare impossible to be NaN")
            });

            if f64::abs(vec[0].0 - vec[1].0) < threshold {
                warn!(
                    "oI lattice, but the two longest vectors {} and {} have almost the same length",
                    vec[0].1, vec[1].1,
                );
            }
            match vec[1].1 {
                Face::A => ExtBravaisClass::oI2,
                Face::B => ExtBravaisClass::oI3,
                Face::C => ExtBravaisClass::oI1,
            }
        }
        BravaisClass::tP => ExtBravaisClass::tP1,
        BravaisClass::tI => {
            if (a - c).abs() < threshold {
                warn!("tI lattice, but a ~ c");
            }

            if c < a {
                ExtBravaisClass::tI1
            } else {
                ExtBravaisClass::tI2
            }
        }
        BravaisClass::hR => {
            if f64::abs(f64::sqrt(3.0) * a - f64::sqrt(2.0) * c) < threshold {
                warn!("hR lattice, but sqrt(3)a almost equal to sqrt(2)c");
            }
            if f64::sqrt(3.0) * a < f64::sqrt(2.0) * c {
                ExtBravaisClass::hR1
            } else {
                ExtBravaisClass::hR2
            }
        }
        BravaisClass::hP => {
            // 143..=163 without 150, 152, 154..=156.
            if [
                143, 144, 145, 146, 147, 148, 149, 151, 153, 157, 159, 160, 161, 162, 163,
            ]
            .contains(&spg_number)
            {
                ExtBravaisClass::hP1
            } else {
                ExtBravaisClass::hP2
            }
        }
        BravaisClass::cP => match spg_number {
            x if (195..=206).contains(&x) => ExtBravaisClass::cP1,
            x if (207..=230).contains(&x) => ExtBravaisClass::cP2,
            _ => unreachable!("cP bravais lattice spacegroup number in wrong range"),
        },
        BravaisClass::cF => match spg_number {
            x if (195..=206).contains(&x) => ExtBravaisClass::cF1,
            x if (207..=230).contains(&x) => ExtBravaisClass::cF2,
            _ => unreachable!("cF bravais lattice spacegroup number in wrong range"),
        },
        BravaisClass::cI => ExtBravaisClass::cI1,
    };

    let path_info = path::lookup(&ext_bravais);
    let path_eval = path::eval(path_info, lattice_params);

    Ok((path_info, crystal_priv))
}
