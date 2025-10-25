use commat::{Atom, CrystalBuilder, Lattice};

fn main() {
    let lattice = Lattice::default();
    let atoms = vec![Atom::new([0.0, 0.0, 0.0], 8)];
    let _ = CrystalBuilder::new()
        .with_lattice(lattice)
        .with_atoms(&atoms)
        .with_atoms(&atoms)
        .build()
        .unwrap();
}
