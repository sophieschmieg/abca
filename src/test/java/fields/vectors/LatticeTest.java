package fields.vectors;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.floatingpoint.FiniteRealVectorSpace;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.integers.Rationals;
import fields.numberfields.NumberField;
import fields.numberfields.NumberFieldIntegers;

class LatticeTest {

	@Test
	void smallTest() {
		Reals r = Reals.r(128);
		FiniteRealVectorSpace space = new FiniteRealVectorSpace(r, 2);
		List<Vector<Real>> basis = new ArrayList<>();
		basis.add(new Vector<>(r.getInteger(4), r.one()));
		basis.add(new Vector<>(r.getInteger(11), r.zero()));
		RealLattice lattice = new RealLattice(space, basis);
		System.out.println(space.latticeReduction(basis, lattice));
	}

	@Test
	void largeTest() {
		int dimension = 20;
		Reals r = Reals.r(128);
		FiniteRealVectorSpace space = new FiniteRealVectorSpace(r, dimension);
		List<Vector<Real>> basis = new ArrayList<>();
		for (int i = 0; i < dimension; i++) {
			basis.add(space.getRandomElement());
		}
		System.out.println(basis);
		RealLattice lattice = new RealLattice(space, basis);
		System.out.println(lattice.getModuleGenerators());
	}

	@Test
	void largeNumberFieldTest() throws IOException {
		Rationals q = Rationals.q();
		NumberField nf = NumberField.getNumberField(q.getUnivariatePolynomialRing().parse("X^32 + 1"));
		NumberFieldIntegers order = nf.maximalOrder();
		System.out.println(order.getModuleGenerators());
		System.out.println(order.idealsOver(2));
		System.out.println(order.idealsOver(3));
		System.out.println(order.idealsOver(5));
		System.out.println(order.idealsOver(7));
		System.out.println(order.idealsOver(11));
	}

}
