package fields.integers;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.integers.Integers.IntE;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Element;
import fields.interfaces.Lattice;
import fields.interfaces.ValueField;
import fields.vectors.AbstractRealInnerProductSpace;
import fields.vectors.EmbeddedRationalLattice;
import fields.vectors.FiniteVectorSpace;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;

public class FiniteRationalVectorSpace extends AbstractRealInnerProductSpace<Fraction, Vector<Fraction>> {
	private Rationals q;
	private FiniteVectorSpace<Fraction> vectorSpace;
	private int dimension;

	public FiniteRationalVectorSpace(int dimension) {
		this.q = Rationals.q();
		this.dimension = dimension;
		this.vectorSpace = new FiniteVectorSpace<>(this.q, this.dimension);
	}

	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	@Override
	public Fraction innerProduct(Vector<Fraction> s1, Vector<Fraction> s2) {
		Fraction result = q.zero();
		for (int i = 0; i < dimension; i++) {
			result = q.add(result, q.multiply(s1.get(i + 1), s2.get(i + 1)));
		}
		return result;
	}

	@Override
	public IntE round(Fraction t) {
		return t.round();
	}

	@Override
	public Fraction fromReal(Real r) {
		return Reals.r(128).roundToFraction(r, 64);
	}

	@Override
	public Real asReal(Fraction t) {
		return Reals.r(128).getFraction(t);
	}

	@Override
	public ValueField<Fraction> getValueField() {
		return q.withInfValue();
	}

	@Override
	public FiniteRationalVectorSpace withDimension(int dimension) {
		return new FiniteRationalVectorSpace(dimension);
	}

	@Override
	public <R extends Element<R>> List<R> latticeReduction(List<R> sublatticeBasis,
			Lattice<R, Fraction, Vector<Fraction>> lattice, double deltaAsDouble) {
		Reals r = getValueField().getReals();
		List<Vector<Fraction>> embeddedLatticeBasis = new ArrayList<>();
		for (R basisVector : sublatticeBasis) {
			embeddedLatticeBasis.add(lattice.embedding(basisVector));
		}
		Matrix<Fraction> originalBaseChange = Matrix.fromColumns(embeddedLatticeBasis);
		MatrixModule<Fraction> matrixModule = originalBaseChange.getModule(q);
		EmbeddedRationalLattice rationalLattice = new EmbeddedRationalLattice(r, dimension(), embeddedLatticeBasis,
				deltaAsDouble);
		List<R> result = new ArrayList<>();
		for (Vector<Fraction> reduced : rationalLattice.getModuleGenerators()) {
			Vector<IntE> asIntegerVector = Vector.mapVector(q.getAsIntegerMap(),
					matrixModule.solve(originalBaseChange, reduced));
			R basisVector = lattice.zero();
			for (int i = 0; i < sublatticeBasis.size(); i++) {
				basisVector = lattice.add(lattice.scalarMultiply(asIntegerVector.get(i + 1), sublatticeBasis.get(i)),
						basisVector);
			}
			result.add(basisVector);
		}
		return result;
	}

	@Override
	public <R extends Element<R>> R closestLatticePoint(Vector<Fraction> t,
			Lattice<R, Fraction, Vector<Fraction>> lattice, double delta) {
		Matrix<Fraction> asMatrix = lattice.generatorsAsMatrix();
		MatrixModule<Fraction> matrixModule = asMatrix.getModule(q);
		Vector<Fraction> asVector = asVector(t);
		int rank = matrixModule.rank(asMatrix);
		if (rank < dimension()) {
			Vector<Fraction> target = asVector;
			asVector = zero();
			for (int i = 0; i < rank; i++) {
				Vector<Fraction> basisVector = asMatrix.column(i + 1);
				Vector<Fraction> projected = scalarMultiply(
						q.divide(innerProduct(target, basisVector), innerProduct(basisVector, basisVector)),
						basisVector);
				asVector = add(projected, asVector);
				target = subtract(target, projected);
			}
		}
		Vector<Fraction> solved = matrixModule.solve(asMatrix, asVector);
		List<IntE> integerSolved = new ArrayList<>();
		for (Fraction coefficient : solved.asList()) {
			integerSolved.add(round(coefficient));
		}
		return lattice.fromVector(new Vector<>(integerSolved));
	}

	@Override
	public Vector<Fraction> zero() {
		return vectorSpace.zero();
	}

	@Override
	public Vector<Fraction> add(Vector<Fraction> s1, Vector<Fraction> s2) {
		return vectorSpace.add(s1, s2);
	}

	@Override
	public Vector<Fraction> negative(Vector<Fraction> s) {
		return vectorSpace.negative(s);
	}

	@Override
	public Vector<Fraction> scalarMultiply(Fraction t, Vector<Fraction> s) {
		return vectorSpace.scalarMultiply(t, s);
	}

	@Override
	public Vector<Fraction> getRandomElement() {
		return vectorSpace.getRandomElement();
	}

	@Override
	public Iterator<Vector<Fraction>> iterator() {
		return vectorSpace.iterator();
	}

	@Override
	public Vector<Fraction> getUnitVector(int index) {
		return vectorSpace.getUnitVector(index);
	}

	@Override
	public List<Vector<Fraction>> getBasis() {
		return vectorSpace.getBasis();
	}

	@Override
	public boolean isLinearIndependent(List<Vector<Fraction>> s) {
		return vectorSpace.isLinearIndependent(s);
	}

	@Override
	public boolean isGeneratingModule(List<Vector<Fraction>> s) {
		return vectorSpace.isGeneratingModule(s);
	}

	@Override
	public List<Vector<Fraction>> getModuleGenerators() {
		return vectorSpace.getModuleGenerators();
	}

	@Override
	public Vector<Fraction> asVector(Vector<Fraction> s) {
		return s;
	}

//	@Override
//	public List<Vector<Fraction>> latticeReduction(List<Vector<Fraction>> s) {
//		return latticeReduction(this, new MathMap<>() {
//			@Override
//			public Fraction evaluate(Fraction t) {
//				return q.getInteger(t.round());
//			}}, s);
//	}
//
//	public List<Vector<IntE>> integerLatticeReduction(List<Vector<IntE>> s) {
//		List<Vector<Fraction>> asFraction = new ArrayList<>();
//		for (Vector<IntE> v : s) {
//			List<Fraction> c = new ArrayList<>();
//			for (IntE value : v.asList()) {
//				c.add(q.getEmbedding(value));
//			}
//			asFraction.add(new Vector<>(c));
//		}
//		asFraction = latticeReduction(asFraction);
//		List<Vector<IntE>> result = new ArrayList<>();
//		for (Vector<Fraction> v : asFraction) {
//			List<IntE> c = new ArrayList<>();
//			for (Fraction value : v.asList()) {
//				c.add(value.asInteger());
//			}
//			result.add(new Vector<>(c));
//		}
//		return result;
//	}

	@Override
	public MatrixAlgebra<Fraction> matrixAlgebra() {
		return vectorSpace.matrixAlgebra();
	}

	@Override
	public int dimension() {
		return dimension;
	}

	@Override
	public List<Vector<Fraction>> nonTrivialCombinations(List<Vector<Fraction>> s) {
		return vectorSpace.nonTrivialCombinations(s);
	}

}
