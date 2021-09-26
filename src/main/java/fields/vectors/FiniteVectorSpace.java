package fields.vectors;

import java.math.BigInteger;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.AbstractModule;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Ring;
import fields.interfaces.VectorSpace;

public class FiniteVectorSpace<T extends Element<T>> extends AbstractModule<T, Vector<T>>
		implements VectorSpace<T, Vector<T>> {
	private FreeModule<T> asModule;
	private Field<T> field;

	public FiniteVectorSpace(Field<T> field, int dimension) {
		this.field = field;
		this.asModule = new FreeModule<>(field, dimension);
	}

	@Override
	public Exactness exactness() {
		return field.exactness();
	}

	@Override
	public Ring<T> getRing() {
		return field;
	}

	@Override
	public Vector<T> zero() {
		return asModule.zero();
	}

	@Override
	public Vector<T> add(Vector<T> s1, Vector<T> s2) {
		return asModule.add(s1, s2);
	}

	@Override
	public Vector<T> negative(Vector<T> s) {
		return asModule.negative(s);
	}

	@Override
	public Vector<T> scalarMultiply(T t, Vector<T> s) {
		return asModule.scalarMultiply(t, s);
	}

	@Override
	public boolean isFree() {
		return true;
	}

	public boolean isBasis(List<Vector<T>> s) {
		return asModule.isBasis(s);
	}
	
	@Override
	public boolean isLinearIndependent(List<Vector<T>> s) {
		return asModule.isLinearIndependent(s);
	}
	
	@Override
	public List<List<T>> nonTrivialCombinations(List<Vector<T>> s) {
		return asModule.nonTrivialCombinations(s);
	}

	@Override
	public boolean isGeneratingModule(List<Vector<T>> s) {
		return asModule.isGeneratingModule(s);
	}

	@Override
	public List<Vector<T>> getModuleGenerators() {
		return asModule.getModuleGenerators();
	}
	
	public Vector<T> getUnitVector(int i) {
		return asModule.getUnitVector(i);
	}

	@Override
	public Vector<T> getRandomElement() {
		return asModule.getRandomElement();
	}

	@Override
	public boolean isFinite() {
		return field.isFinite();
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return asModule.getNumberOfElements();
	}

	@Override
	public Iterator<Vector<T>> iterator() {
		return asModule.iterator();
	}

	@Override
	public Field<T> getField() {
		return field;
	}

	@Override
	public List<Vector<T>> getBasis() {
		return asModule.getBasis();
	}

	@Override
	public Vector<T> asVector(Vector<T> s) {
		return s;
	}
	
	@Override
	public MatrixAlgebra<T> matrixAlgebra() {
		return asModule.matrixAlgebra();
	}
	
	@Override
	public int dimension() {
		return asModule.dimension();
	}
	
	@Override
	public String toString() {
		return field.toString() + "^" + dimension();
	}
}
