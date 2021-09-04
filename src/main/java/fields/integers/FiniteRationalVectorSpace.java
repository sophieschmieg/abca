package fields.integers;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import fields.helper.AbstractInnerProductSpace;
import fields.integers.Integers.IntE;
import fields.integers.Rationals.Fraction;
import fields.interfaces.MathMap;
import fields.interfaces.ValueField;
import fields.vectors.FiniteVectorSpace;
import fields.vectors.MatrixAlgebra;
import fields.vectors.Vector;

public class FiniteRationalVectorSpace extends AbstractInnerProductSpace<Fraction, Vector<Fraction>> {
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
		return Exactness.FLOATING_POINT;
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
	public ValueField<Fraction> getValueField() {
		return q.withInfValue();
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
	
	@Override
	public List<Vector<Fraction>> latticeReduction(List<Vector<Fraction>> s) {
		return latticeReduction(this, new MathMap<>() {
			@Override
			public Fraction evaluate(Fraction t) {
				return q.getInteger(t.round());
			}}, s);
	}

	public List<Vector<IntE>> integerLatticeReduction(List<Vector<IntE>> s) {
		List<Vector<Fraction>> asFraction = new ArrayList<>();
		for (Vector<IntE> v : s) {
			List<Fraction> c = new ArrayList<>();
			for (IntE value : v.asList()) {
				c.add(q.getEmbedding(value));
			}
			asFraction.add(new Vector<>(c));
		}
		asFraction = latticeReduction(asFraction);
		List<Vector<IntE>> result = new ArrayList<>();
		for (Vector<Fraction> v : asFraction) {
			List<IntE> c = new ArrayList<>();
			for (Fraction value : v.asList()) {
				c.add(value.asInteger());
			}
			result.add(new Vector<>(c));
		}
		return result;
	}

	@Override
	public MatrixAlgebra<Fraction> matrixAlgebra() {
		return vectorSpace.matrixAlgebra();
	}
	
	@Override
	public int dimension() {
		return dimension;
	}

	@Override
	public List<List<Fraction>> nonTrivialCombinations(List<Vector<Fraction>> s) {
		return vectorSpace.nonTrivialCombinations(s);
	}

}
