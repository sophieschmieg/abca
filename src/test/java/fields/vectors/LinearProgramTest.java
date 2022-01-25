package fields.vectors;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.floatingpoint.FiniteRealVectorSpace;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.integers.FiniteRationalVectorSpace;
import fields.integers.Integers;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Lattice;
import fields.interfaces.RealInnerProductSpace.LinearProgramResult;
import fields.interfaces.RealInnerProductSpace.SimplexAlgorithmStatus;
import fields.vectors.DualVectorSpace.Dual;
import graphics.Canvas2D;

class LinearProgramTest {

	@Test
	void testRectangle() {
		Rationals q = Rationals.q();
		FiniteRationalVectorSpace space = new FiniteRationalVectorSpace(2);
		DualVectorSpace<Fraction, Vector<Fraction>> dualSpace = space.getDual();
		List<Dual<Fraction, Vector<Fraction>>> constraints = new ArrayList<>();
		List<Fraction> rhs = new ArrayList<>();
		constraints.add(dualSpace.negative(dualSpace.getUnitVector(1)));
		constraints.add(dualSpace.negative(dualSpace.getUnitVector(2)));
		constraints.add(dualSpace.getUnitVector(1));
		constraints.add(dualSpace.getUnitVector(2));
		rhs.add(q.getInteger(-2));
		rhs.add(q.getInteger(-3));
		rhs.add(q.getInteger(4));
		rhs.add(q.getInteger(5));
		List<Fraction> maximize = new ArrayList<>();
		maximize.add(q.getInteger(1));
		maximize.add(q.getInteger(1));
		LinearProgramResult<Fraction, Vector<Fraction>> result = space.linearProgram(
				new Polytope<>(space, constraints, rhs), dualSpace.canonicalIsomorphism(new Vector<>(maximize)));
		assertEquals(SimplexAlgorithmStatus.OKAY, result.getStatus());
		assertEquals(new Vector<>(q.getInteger(4), q.getInteger(5)), result.getSolution());
		assertEquals(q.getInteger(9), result.getValue());
	}

	@Test
	void testOctagon() {
		Rationals q = Rationals.q();
		FiniteRationalVectorSpace space = new FiniteRationalVectorSpace(2);
		DualVectorSpace<Fraction, Vector<Fraction>> dualSpace = space.getDual();
		List<Dual<Fraction, Vector<Fraction>>> constraints = new ArrayList<>();
		List<Fraction> rhs = new ArrayList<>();
		constraints.add(dualSpace.negative(dualSpace.getUnitVector(1)));
		constraints.add(dualSpace.add(dualSpace.negative(dualSpace.getUnitVector(1)),
				dualSpace.negative(dualSpace.getUnitVector(2))));
		constraints.add(dualSpace.negative(dualSpace.getUnitVector(2)));
		constraints.add(dualSpace.add(dualSpace.getUnitVector(1), dualSpace.negative(dualSpace.getUnitVector(2))));
		constraints.add(dualSpace.getUnitVector(1));
		constraints.add(dualSpace.add(dualSpace.getUnitVector(1), dualSpace.getUnitVector(2)));
		constraints.add(dualSpace.getUnitVector(2));
		constraints.add(dualSpace.add(dualSpace.negative(dualSpace.getUnitVector(1)), dualSpace.getUnitVector(2)));
		rhs.add(q.getInteger(2));
		rhs.add(q.getInteger(3));
		rhs.add(q.getInteger(2));
		rhs.add(q.getInteger(3));
		rhs.add(q.getInteger(2));
		rhs.add(q.getInteger(3));
		rhs.add(q.getInteger(2));
		rhs.add(q.getInteger(3));
		List<Fraction> maximize = new ArrayList<>();
		maximize.add(q.getInteger(2));
		maximize.add(q.getInteger(1));
		LinearProgramResult<Fraction, Vector<Fraction>> result = space.linearProgram(
				new Polytope<>(space, constraints, rhs), dualSpace.canonicalIsomorphism(new Vector<>(maximize)));
		assertEquals(SimplexAlgorithmStatus.OKAY, result.getStatus());
		assertEquals(new Vector<>(q.getInteger(2), q.getInteger(1)), result.getSolution());
		assertEquals(q.getInteger(5), result.getValue());
	}

	@Test
	void testOctagonNotFeasibleStart() {
		Rationals q = Rationals.q();
		FiniteRationalVectorSpace space = new FiniteRationalVectorSpace(2);
		DualVectorSpace<Fraction, Vector<Fraction>> dualSpace = space.getDual();
		List<Dual<Fraction, Vector<Fraction>>> constraints = new ArrayList<>();
		List<Fraction> rhs = new ArrayList<>();
		constraints.add(dualSpace.negative(dualSpace.getUnitVector(1)));
		constraints.add(dualSpace.getUnitVector(1));
		constraints.add(dualSpace.getUnitVector(2));
		constraints.add(dualSpace.add(dualSpace.negative(dualSpace.getUnitVector(1)),
				dualSpace.negative(dualSpace.getUnitVector(2))));
		constraints.add(dualSpace.negative(dualSpace.getUnitVector(2)));
		constraints.add(dualSpace.add(dualSpace.getUnitVector(1), dualSpace.negative(dualSpace.getUnitVector(2))));
		constraints.add(dualSpace.add(dualSpace.getUnitVector(1), dualSpace.getUnitVector(2)));
		constraints.add(dualSpace.add(dualSpace.negative(dualSpace.getUnitVector(1)), dualSpace.getUnitVector(2)));
		rhs.add(q.getInteger(2));
		rhs.add(q.getInteger(2));
		rhs.add(q.getInteger(2));
		rhs.add(q.getInteger(3));
		rhs.add(q.getInteger(2));
		rhs.add(q.getInteger(3));
		rhs.add(q.getInteger(3));
		rhs.add(q.getInteger(3));
		List<Fraction> maximize = new ArrayList<>();
		maximize.add(q.getInteger(2));
		maximize.add(q.getInteger(1));
		LinearProgramResult<Fraction, Vector<Fraction>> result = space.linearProgram(
				new Polytope<>(space, constraints, rhs), dualSpace.canonicalIsomorphism(new Vector<>(maximize)));
		assertEquals(SimplexAlgorithmStatus.OKAY, result.getStatus());
		assertEquals(new Vector<>(q.getInteger(2), q.getInteger(1)), result.getSolution());
		assertEquals(q.getInteger(5), result.getValue());
	}

	@Test
	void testIntegerProgram() {
		Reals r = Reals.r(128);
		FiniteRealVectorSpace space = new FiniteRealVectorSpace(r, 2);
		List<Vector<Real>> latticeGenerators = new ArrayList<>();
		latticeGenerators.add(space.getUnitVector(1));
		latticeGenerators.add(
				new Vector<>(r.inverse(r.getInteger(2)), r.divide(r.positiveSqrt(r.getInteger(3)), r.getInteger(2))));
		Lattice<Vector<Real>, Real, Vector<Real>> lattice = new RealLattice(space, latticeGenerators);
		DualVectorSpace<Real, Vector<Real>> dualSpace = space.getDual();
		List<Dual<Real, Vector<Real>>> constraints = new ArrayList<>();
		List<Real> rhs = new ArrayList<>();
		constraints.add(dualSpace.negative(dualSpace.getUnitVector(1)));
		constraints.add(dualSpace.getUnitVector(1));
		constraints.add(dualSpace.getUnitVector(2));
		constraints.add(dualSpace.add(dualSpace.negative(dualSpace.getUnitVector(1)),
				dualSpace.negative(dualSpace.getUnitVector(2))));
		constraints.add(dualSpace.negative(dualSpace.getUnitVector(2)));
		constraints.add(dualSpace.add(dualSpace.getUnitVector(1), dualSpace.negative(dualSpace.getUnitVector(2))));
		constraints.add(dualSpace.add(dualSpace.getUnitVector(1), dualSpace.getUnitVector(2)));
		constraints.add(dualSpace.add(dualSpace.negative(dualSpace.getUnitVector(1)), dualSpace.getUnitVector(2)));
		rhs.add(r.getInteger(2));
		rhs.add(r.getInteger(2));
		rhs.add(r.getInteger(2));
		rhs.add(r.getInteger(3));
		rhs.add(r.getInteger(2));
		rhs.add(r.getInteger(3));
		rhs.add(r.getInteger(3));
		rhs.add(r.getInteger(3));
		List<Real> maximize = new ArrayList<>();
		maximize.add(r.getInteger(1));
		maximize.add(r.getInteger(1));
		Vector<Real> result = space.latticeLinearProgram(new Polytope<>(space, constraints, rhs),
				dualSpace.canonicalIsomorphism(new Vector<>(maximize)), lattice);
		System.out.println(result);
		System.out.println(lattice.asVector(result));
		Integers z = Integers.z();
		assertEquals(new Vector<>(z.zero(), z.getInteger(2)), lattice.asVector(result));
		List<Vector<Real>> latticeVerteces = space
				.latticeVertexPointsInPolytope(new Polytope<>(space, constraints, rhs), lattice);
		System.out.println(latticeVerteces);
		List<Vector<Real>> latticePoints = space.latticePointsInPolytope(new Polytope<>(space, constraints, rhs),
				lattice);
		System.out.println(latticePoints);
	}

	@Test
	void testIntegerPolygon() {
		Reals r = Reals.r(128);
		FiniteRealVectorSpace space = new FiniteRealVectorSpace(r, 2);
		List<Vector<Real>> latticeGenerators = new ArrayList<>();
		latticeGenerators.add(space.getUnitVector(1));
		latticeGenerators.add(
				new Vector<>(r.inverse(r.getInteger(2)), r.divide(r.positiveSqrt(r.getInteger(3)), r.getInteger(2))));
		Lattice<Vector<Real>, Real, Vector<Real>> lattice = new RealLattice(space, latticeGenerators);
		DualVectorSpace<Real, Vector<Real>> dualSpace = space.getDual();
		List<Dual<Real, Vector<Real>>> constraints = new ArrayList<>();
		List<Real> rhs = new ArrayList<>();
		constraints.add(dualSpace.negative(dualSpace.getUnitVector(1)));
		constraints.add(dualSpace.getUnitVector(1));
		constraints.add(dualSpace.getUnitVector(2));
		constraints.add(dualSpace.add(dualSpace.negative(dualSpace.getUnitVector(1)),
				dualSpace.negative(dualSpace.getUnitVector(2))));
		constraints.add(dualSpace.negative(dualSpace.getUnitVector(2)));
		constraints.add(dualSpace.add(dualSpace.getUnitVector(1), dualSpace.negative(dualSpace.getUnitVector(2))));
		constraints.add(dualSpace.add(dualSpace.getUnitVector(1), dualSpace.getUnitVector(2)));
		constraints.add(dualSpace.add(dualSpace.negative(dualSpace.getUnitVector(1)), dualSpace.getUnitVector(2)));
		rhs.add(r.getInteger(2));
		rhs.add(r.getInteger(2));
		rhs.add(r.getInteger(2));
		rhs.add(r.getInteger(3));
		rhs.add(r.getInteger(2));
		rhs.add(r.getInteger(3));
		rhs.add(r.getInteger(3));
		rhs.add(r.getInteger(3));
		List<Real> maximize = new ArrayList<>();
		maximize.add(r.getInteger(1));
		maximize.add(r.getInteger(1));
		Polytope<Real, Vector<Real>> polytope = new Polytope<>(space, constraints, rhs);
		Canvas2D<Real, Vector<Real>> canvas = new Canvas2D<>(space, 1024, 1024, r.getInteger(-3), r.getInteger(3),
				r.getInteger(-3), r.getInteger(3));
		canvas.setLattice(lattice, -5, 5, -5, 5);
		canvas.setPolytope(polytope);
		try {
			Thread.sleep(10000);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		canvas.setPolytope(Polytope.fromVertices(space, space.latticePointsInPolytope(polytope, lattice)));
		try {
			Thread.sleep(10000);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		//		Vector<Real> solution = space.latticeLinearProgram(polytope, constraints.get(2), lattice);
		List<Vector<Real>> result = space.latticePointsInPolytope(polytope, lattice);
		System.out.println(result);
	}

	@Test
	void testConvexHull() {
		Rationals q = Rationals.q();
		FiniteRationalVectorSpace space = new FiniteRationalVectorSpace(3);
		List<Vector<Fraction>> vertices = new ArrayList<>();
		vertices.add(new Vector<>(q.getInteger(-1), q.getInteger(-2), q.getInteger(1)));
		vertices.add(new Vector<>(q.getInteger(-2), q.getInteger(-1), q.getInteger(0)));
		vertices.add(new Vector<>(q.getInteger(-1), q.getInteger(2), q.getInteger(1)));
		vertices.add(new Vector<>(q.getInteger(-2), q.getInteger(1), q.getInteger(0)));
		vertices.add(new Vector<>(q.getInteger(1), q.getInteger(-2), q.getInteger(3)));
		vertices.add(new Vector<>(q.getInteger(2), q.getInteger(-1), q.getInteger(4)));
		vertices.add(new Vector<>(q.getInteger(1), q.getInteger(2), q.getInteger(3)));
		vertices.add(new Vector<>(q.getInteger(2), q.getInteger(1), q.getInteger(4)));
		Polytope<Fraction, Vector<Fraction>> polytope = Polytope.fromVertices(space, vertices);
		System.out.println(polytope.vertices());
	}

	@Test
	void testVertices() {
		Rationals q = Rationals.q();
		FiniteRationalVectorSpace space = new FiniteRationalVectorSpace(2);
		DualVectorSpace<Fraction, Vector<Fraction>> dualSpace = space.getDual();
		List<Dual<Fraction, Vector<Fraction>>> constraints = new ArrayList<>();
		List<Fraction> rhs = new ArrayList<>();
		constraints.add(dualSpace.negative(dualSpace.getUnitVector(1)));
		constraints.add(dualSpace.getUnitVector(1));
		constraints.add(dualSpace.getUnitVector(2));
		constraints.add(dualSpace.add(dualSpace.negative(dualSpace.getUnitVector(1)),
				dualSpace.negative(dualSpace.getUnitVector(2))));
		constraints.add(dualSpace.negative(dualSpace.getUnitVector(2)));
		constraints.add(dualSpace.add(dualSpace.getUnitVector(1), dualSpace.negative(dualSpace.getUnitVector(2))));
		constraints.add(dualSpace.add(dualSpace.getUnitVector(1), dualSpace.getUnitVector(2)));
		constraints.add(dualSpace.add(dualSpace.negative(dualSpace.getUnitVector(1)), dualSpace.getUnitVector(2)));
		rhs.add(q.getInteger(2));
		rhs.add(q.getInteger(2));
		rhs.add(q.getInteger(2));
		rhs.add(q.getInteger(3));
		rhs.add(q.getInteger(2));
		rhs.add(q.getInteger(3));
		rhs.add(q.getInteger(3));
		rhs.add(q.getInteger(3));
		Polytope<Fraction, Vector<Fraction>> polytope = new Polytope<>(space, constraints, rhs);
		System.out.println(polytope.vertices());
	}

	@Test
	void testOctahedron() {
		Rationals q = Rationals.q();
		FiniteRationalVectorSpace space = new FiniteRationalVectorSpace(3);
		List<Vector<Fraction>> vertices = new ArrayList<>();
		vertices.add(new Vector<>(q.getInteger(-1), q.getInteger(-1), q.getInteger(0)));
		vertices.add(new Vector<>(q.getInteger(-1), q.getInteger(1), q.getInteger(0)));
		vertices.add(new Vector<>(q.getInteger(1), q.getInteger(-1), q.getInteger(0)));
		vertices.add(new Vector<>(q.getInteger(1), q.getInteger(1), q.getInteger(0)));
		vertices.add(new Vector<>(q.getInteger(0), q.getInteger(0), q.getInteger(1)));
		vertices.add(new Vector<>(q.getInteger(0), q.getInteger(0), q.getInteger(-1)));
		Polytope<Fraction, Vector<Fraction>> polytope = Polytope.fromVertices(space, vertices);
		Polytope<Fraction, Vector<Fraction>> polytopeCopy = new Polytope<>(space, polytope.getConstraints(),
				polytope.getRightHandSide());
		System.out.println(polytopeCopy.vertices());

	}
}
