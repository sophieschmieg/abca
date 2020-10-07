package fields.floatingpoint;

import java.util.Iterator;
import java.util.List;

import fields.floatingpoint.Complex.ComplexNumber;
import fields.helper.AbstractInnerProductSpace;
import fields.interfaces.ValueField;
import fields.vectors.FreeModule;
import fields.vectors.Vector;

public class FiniteComplexVectorSpace extends AbstractInnerProductSpace<ComplexNumber, Vector<ComplexNumber>> {
	private Complex c;
	private FreeModule<ComplexNumber> freeModule;
	private int dimension;

	public FiniteComplexVectorSpace(int dimension) {
		this.c = Complex.c();
		this.dimension = dimension;
		this.freeModule = new FreeModule<>(this.c, this.dimension);
	}

	@Override
	public ComplexNumber innerProduct(Vector<ComplexNumber> s1, Vector<ComplexNumber> s2) {
		ComplexNumber result = c.zero();
		for (int i = 0; i < dimension; i++) {
			result = c.add(result, c.multiply(s1.get(i + 1), c.conjugate(s2.get(i + 1))));
		}
		return result;
	}

	@Override
	public ValueField<ComplexNumber> getValueField() {
		return c;
	}

	@Override
	public Vector<ComplexNumber> zero() {
		return freeModule.zero();
	}

	@Override
	public Vector<ComplexNumber> add(Vector<ComplexNumber> s1, Vector<ComplexNumber> s2) {
		return freeModule.add(s1, s2);
	}

	@Override
	public Vector<ComplexNumber> negative(Vector<ComplexNumber> s) {
		return freeModule.negative(s);
	}

	@Override
	public Vector<ComplexNumber> scalarMultiply(ComplexNumber t, Vector<ComplexNumber> s) {
		return freeModule.scalarMultiply(t, s);
	}

	@Override
	public Vector<ComplexNumber> getRandomElement() {
		return freeModule.getRandomElement();
	}

	@Override
	public Iterator<Vector<ComplexNumber>> iterator() {
		return freeModule.iterator();
	}

	@Override
	public List<Vector<ComplexNumber>> getBasis() {
		return freeModule.getBasis();
	}

	@Override
	public boolean isLinearIndependent(List<Vector<ComplexNumber>> s) {
		return freeModule.isLinearIndependent(s);
	}

	@Override
	public boolean isGeneratingModule(List<Vector<ComplexNumber>> s) {
		return freeModule.isGeneratingModule(s);
	}

	@Override
	public List<Vector<ComplexNumber>> getModuleGenerators() {
		return freeModule.getModuleGenerators();
	}

	@Override
	public Vector<ComplexNumber> asVector(Vector<ComplexNumber> s) {
		return s;
	}
}
