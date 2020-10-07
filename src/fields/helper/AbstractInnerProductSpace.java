package fields.helper;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.InnerProductSpace;
import fields.interfaces.Ring;
import fields.interfaces.ValueField;

public abstract class AbstractInnerProductSpace<T extends Element<T>, S extends Element<S>> extends AbstractModule<T, S>
		implements InnerProductSpace<T, S> {

	@Override
	public boolean isFree() {
		return true;
	}

	@Override
	public boolean isFinite() {
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		throw new InfinityException();
	}

	public double norm(S s) {
		return Math.sqrt(getValueField().value(innerProduct(s, s)));
	}

	public Ring<T> getRing() {
		return getValueField();
	}

	public Field<T> getField() {
		return getValueField();
	}

	public List<S> getGenerators() {
		return getBasis();
	}

	@Override
	public List<S> gramSchmidt(List<S> s) {
		List<S> orthogonal = new ArrayList<>();
		ValueField<T> f = getValueField();
		for (int i = 0; i < s.size(); i++) {
			S vector = s.get(i);
			for (int j = 0; j < i; j++) {
				vector = subtract(vector, scalarMultiply(
						f.divide(innerProduct(s.get(i), s.get(j)), innerProduct(s.get(j), s.get(j))), s.get(j)));
			}
			orthogonal.add(vector);
		}
		return orthogonal;
	}

	@Override
	public List<S> normedGramSchmidt(List<S> s) {
		List<S> orthogonal = gramSchmidt(s);
		List<S> orthonormal = new ArrayList<>();
		ValueField<T> f = getValueField();
		for (S vector : orthogonal) {
			orthonormal.add(scalarMultiply(f.inverse(f.sqrt(innerProduct(vector, vector)).get(0)), vector));
		}
		return orthonormal;
	}

}
