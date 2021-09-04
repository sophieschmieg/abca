package fields.floatingpoint;

import java.util.Iterator;
import java.util.List;

import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractInnerProductSpace;
import fields.interfaces.MathMap;
import fields.interfaces.ValueField;
import fields.vectors.FiniteVectorSpace;
import fields.vectors.MatrixAlgebra;
import fields.vectors.Vector;

public class FiniteRealVectorSpace extends AbstractInnerProductSpace<Real, Vector<Real>> {
	private Reals r;
	private FiniteVectorSpace<Real> vectorSpace;
	private int dimension;

	public FiniteRealVectorSpace(Reals r, int dimension) {
		this.r = r;
		this.dimension = dimension;
		this.vectorSpace = new FiniteVectorSpace<>(this.r, this.dimension);
	}

	@Override
	public Exactness exactness() {
		return Exactness.FLOATING_POINT;
	}

	@Override
	public Real innerProduct(Vector<Real> s1, Vector<Real> s2) {
		Real result = r.zero();
		for (int i = 0; i < dimension; i++) {
			result = r.add(result, r.multiply(s1.get(i + 1), s2.get(i + 1)));
		}
		return result;
	}

	@Override
	public ValueField<Real> getValueField() {
		return r;
	}

	@Override
	public Vector<Real> zero() {
		return vectorSpace.zero();
	}

	@Override
	public Vector<Real> add(Vector<Real> s1, Vector<Real> s2) {
		return vectorSpace.add(s1, s2);
	}

	@Override
	public Vector<Real> negative(Vector<Real> s) {
		return vectorSpace.negative(s);
	}

	@Override
	public Vector<Real> scalarMultiply(Real t, Vector<Real> s) {
		return vectorSpace.scalarMultiply(t, s);
	}

	@Override
	public Vector<Real> getRandomElement() {
		return vectorSpace.getRandomElement();
	}

	@Override
	public Iterator<Vector<Real>> iterator() {
		return vectorSpace.iterator();
	}

	@Override
	public List<Vector<Real>> getBasis() {
		return vectorSpace.getBasis();
	}

	@Override
	public boolean isLinearIndependent(List<Vector<Real>> s) {
		return vectorSpace.isLinearIndependent(s);
	}

	@Override
	public boolean isGeneratingModule(List<Vector<Real>> s) {
		return vectorSpace.isGeneratingModule(s);
	}

	@Override
	public List<Vector<Real>> getModuleGenerators() {
		return vectorSpace.getModuleGenerators();
	}

	@Override
	public Vector<Real> asVector(Vector<Real> s) {
		return s;
	}
	
	@Override
	public List<Vector<Real>> latticeReduction(List<Vector<Real>> s) {
		return latticeReduction(this, new MathMap<>() {
			@Override
			public Real evaluate(Real t) {
				return r.getInteger(t.round());
			}}, s);
	}

	@Override
	public MatrixAlgebra<Real> matrixAlgebra() {
		return vectorSpace.matrixAlgebra();
	}
	
	@Override
	public int dimension() {
		return dimension;
	}

	@Override
	public List<List<Real>> nonTrivialCombinations(List<Vector<Real>> s) {
		return vectorSpace.nonTrivialCombinations(s);
	}

}
