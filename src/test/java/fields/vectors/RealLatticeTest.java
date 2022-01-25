package fields.vectors;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.floatingpoint.FiniteRealVectorSpace;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;

class RealLatticeTest {

	@Test
	void test2() {
		Reals r = Reals.r(128);
		FiniteRealVectorSpace space = new FiniteRealVectorSpace(r, 2);
		List<Vector<Real>> generators = new ArrayList<>();
		generators.add(new Vector<>(r.getInteger(6), r.getInteger(4)));
		generators.add(new Vector<>(r.getInteger(9), r.getInteger(6)));
		RealLattice lattice = new RealLattice(space, generators);
		assertEquals(1, lattice.rank());
		for (Vector<Real> vector : generators) {
			assertTrue(lattice.contains(vector));
		}
	}

	@Test
	void test3() {
		Reals r = Reals.r(128);
		FiniteRealVectorSpace space = new FiniteRealVectorSpace(r, 2);
		List<Vector<Real>> generators = new ArrayList<>();
		generators.add(new Vector<>(r.getInteger(5), r.getInteger(10)));
		generators.add(new Vector<>(r.getInteger(0), r.getInteger(-5)));
		generators.add(new Vector<>(r.getInteger(-3), r.getInteger(2)));
		RealLattice lattice = new RealLattice(space, generators);
		assertEquals(2, lattice.rank());
		for (Vector<Real> vector : generators) {
			assertTrue(lattice.contains(vector));
		}
	}

	@Test
	void test31() {
		Reals r = Reals.r(128);
		FiniteRealVectorSpace space = new FiniteRealVectorSpace(r, 2);
		List<Vector<Real>> generators = new ArrayList<>();
		generators.add(new Vector<>(r.getInteger(6), r.getInteger(-4)));
		generators.add(new Vector<>(r.getDouble(4.5), r.getInteger(-3)));
		generators.add(new Vector<>(r.getInteger(-3), r.getInteger(2)));
		RealLattice lattice = new RealLattice(space, generators);
		assertEquals(1, lattice.rank());
		for (Vector<Real> vector : generators) {
			assertTrue(lattice.contains(vector));
		}
	}

	@Test
	void testPi() {
		Reals r = Reals.r(128);
		FiniteRealVectorSpace space = new FiniteRealVectorSpace(r, 2);
		List<Vector<Real>> generators = new ArrayList<>();
		generators.add(new Vector<>(r.pi(), r.getInteger(10)));
		generators.add(new Vector<>(r.getInteger(-3), r.getInteger(2)));
		RealLattice lattice = new RealLattice(space, generators);
		assertEquals(2, lattice.rank());
		for (Vector<Real> vector : generators) {
			assertTrue(lattice.contains(vector));
		}
	}

}
