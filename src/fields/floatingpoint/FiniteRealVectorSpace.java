package fields.floatingpoint;

import java.util.Iterator;
import java.util.List;

import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractInnerProductSpace;
import fields.interfaces.ValueField;
import fields.vectors.FreeModule;
import fields.vectors.Vector;

public class FiniteRealVectorSpace extends AbstractInnerProductSpace<Real, Vector<Real>> {
	private Reals r;
	private FreeModule<Real> freeModule;
	private int dimension;

	public FiniteRealVectorSpace(int dimension) {
		this.r = Reals.r();
		this.dimension = dimension;
		this.freeModule = new FreeModule<>(this.r, this.dimension);
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
		return freeModule.zero();
	}

	@Override
	public Vector<Real> add(Vector<Real> s1, Vector<Real> s2) {
		return freeModule.add(s1, s2);
	}

	@Override
	public Vector<Real> negative(Vector<Real> s) {
		return freeModule.negative(s);
	}

	@Override
	public Vector<Real> scalarMultiply(Real t, Vector<Real> s) {
		return freeModule.scalarMultiply(t, s);
	}

	@Override
	public Vector<Real> getRandomElement() {
		return freeModule.getRandomElement();
	}

	@Override
	public Iterator<Vector<Real>> iterator() {
		return freeModule.iterator();
	}

	@Override
	public List<Vector<Real>> getBasis() {
		return freeModule.getBasis();
	}

	@Override
	public boolean isLinearIndependent(List<Vector<Real>> s) {
		return freeModule.isLinearIndependent(s);
	}

	@Override
	public boolean isGeneratingModule(List<Vector<Real>> s) {
		return freeModule.isGeneratingModule(s);
	}

	@Override
	public List<Vector<Real>> getModuleGenerators() {
		return freeModule.getModuleGenerators();
	}

	@Override
	public Vector<Real> asVector(Vector<Real> s) {
		return s;
	}
}
