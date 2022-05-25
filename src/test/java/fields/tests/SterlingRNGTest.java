package fields.tests;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.floatingpoint.FiniteRealVectorSpace;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.vectors.RealLattice;
import fields.vectors.SubRealInnerProductSpace;
import fields.vectors.Vector;

class SterlingRNGTest {

	@Test
	void test() {
		int precision = 1024;
		int bitRate = 256;
		Reals r = Reals.r(precision);
		Integers z = Integers.z();
		IntE n = z.getRandomElement(z.power(z.getInteger(2), bitRate));
		Real sqrt = r.positiveSqrt(r.getInteger(n));
		System.out.println("Secret integer: " + n);
		System.out.println("Square root: " + sqrt);
		System.out.println("Fractional part: " + sqrt.fractionalPart());
		IntE reconstructed = breakBadRng(sqrt.fractionalPart(), precision, bitRate);
		System.out.println("Reconstructed integer: " + reconstructed);
		assertEquals(n, reconstructed);
	}

	private IntE breakBadRng(Real sqrt, int precision, int bitRate) {
		Reals r = Reals.r(precision);
		Integers z = Integers.z();
		FiniteRealVectorSpace space = new FiniteRealVectorSpace(r, 3);
		Real power = r.getPowerOfTwo(bitRate + 4);
		List<Vector<Real>> generators = new ArrayList<>();
		generators.add(new Vector<>(r.one(), r.zero(), r.multiply(2, sqrt, power)));
		generators.add(new Vector<>(r.zero(), r.one(), r.negative(power)));
		for (Vector<Real> generator : generators) {
			System.out.println("Lattice generator: " + generator);
		}
		SubRealInnerProductSpace<Real, Vector<Real>> subSpace = new SubRealInnerProductSpace<>(space, generators);
		RealLattice lattice = new RealLattice(space, generators);
		System.out.println("LLL reducing lattice basis...");
		for (Vector<Real> generator : lattice.getModuleGenerators()) {
			System.out.println("Lattice generator: " + generator);
		}
		Vector<Real> target = new Vector<>(r.zero(), r.zero(), r.multiply(-1, power, sqrt, sqrt));
		System.out.println("Target vector: " + target);
		target = subSpace.project(target);
		System.out.println("Target projected to lattice span: " + target);
		Vector<IntE> result = lattice.asVector(space.closestLatticePoint(target, lattice));
		Vector<Real> closest = lattice.fromVector(result);
		System.out.println("Closest lattice vector: " + closest);
		IntE intPart = closest.get(1).round();
		System.out.println("Reconstructed integer part: " + intPart);
		IntE diff = closest.get(2).round();
		System.out.println("Reconstructed difference to integer part square: " + diff);
		return z.add(z.multiply(intPart, intPart), diff);
	}
}
