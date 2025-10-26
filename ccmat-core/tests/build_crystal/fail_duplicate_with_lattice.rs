use ccmat::{lattice_angstrom, Atom, CrystalBuilder};

fn main() {
    let lattice = lattice_angstrom![
        a = (1.0, 0.0, 0.0),
        b = (0.0, 1.0, 0.0),
        c = (0.0, 0.0, 1.0),
    ];
    let atoms = vec![Atom::new([0.0, 0.0, 0.0], 8)];
    let _ = CrystalBuilder::new()
        .with_lattice(lattice)
        .with_lattice(lattice)
        .with_atoms(&atoms)
        .build()
        .unwrap();
}
