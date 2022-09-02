package fields.local;

import java.math.BigInteger;
import java.util.Iterator;

import fields.exceptions.InfinityException;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractField;
import fields.interfaces.DiscreteValuationField;
import fields.interfaces.DiscreteValuationRing;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.UnivariatePolynomial;
import util.Identity;

public class TrivialDiscreteValuationField<T extends Element<T>> extends AbstractField<T>
		implements DiscreteValuationField<T, T> {
	private Field<T> field;

	public TrivialDiscreteValuationField(Field<T> field) {
		this.field = field;
	}

	@Override
	public String toString() {
		return field.toString();
	}

	@Override
	public Real value(T t) {
		if (t.equals(field.zero())) {
			return Reals.r(1024).zero();
		}
		return Reals.r(1024).one();
	}
	
	@Override
	public Reals getReals() {
		return Reals.r(1024);
	}

	@Override
	public T zero() {
		return field.zero();
	}

	@Override
	public T one() {
		return field.one();
	}

	@Override
	public BigInteger characteristic() {
		return field.characteristic();
	}

	@Override
	public T add(T t1, T t2) {
		return field.add(t1, t2);
	}

	@Override
	public T negative(T t) {
		return field.negative(t);
	}

	@Override
	public T multiply(T t1, T t2) {
		return field.multiply(t1, t2);
	}

	@Override
	public T inverse(T t) {
		return field.inverse(t);
	}

	@Override
	public Exactness exactness() {
		return field.exactness();
	}

	@Override
	public T getRandomElement() {
		return field.getRandomElement();
	}

	@Override
	public boolean isFinite() {
		return field.isFinite();
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return field.getNumberOfElements();
	}

	@Override
	public Iterator<T> iterator() {
		return field.iterator();
	}

	@Override
	public boolean isComplete() {
		return true;
	}

	@Override
	public T inverse(T t, int accuracy) {
		return inverse(t);
	}

	@Override
	public T negative(T t, int accuracy) {
		return negative(t);
	}

	@Override
	public Value valuation(T t) {
		if (t.equals(field.zero())) {
			return Value.INFINITY;
		}
		return Value.ZERO;
	}

	@Override
	public T uniformizer() {
		return one();
	}

	@Override
	public Field<T> residueField() {
		return this;
	}

	@Override
	public T reduceInteger(T t) {
		return t;
	}

	@Override
	public T upToUniformizerPower(T t) {
		return t;
	}

	@Override
	public T liftToInteger(T s) {
		return s;
	}

	@Override
	public T round(T t, int accuracy) {
		return t;
	}

	@Override
	public int getAccuracy() {
		return 1;
	}

	@Override
	public DiscreteValuationField<T, T> withAccuracy(int accuracy) {
		return this;
	}

	@Override
	public OtherVersion<T, T, T, TrivialDiscreteValuationField<T>> exact() {
		return new OtherVersion<>(this, new Identity<>(), new Identity<>());
	}

	@Override
	public OtherVersion<T, T, T, TrivialDiscreteValuationField<T>> complete(int accuracy) {
		return new OtherVersion<>(this, new Identity<>(), new Identity<>());
	}

	@Override
	public Extension<T, ?, ?, ?> getExtension(UnivariatePolynomial<T> minimalPolynomial) {
		return field.getExtension(minimalPolynomial);
	}

	@Override
	public ExtensionOfDiscreteValuationField<T, T, ?, ?, ?, ?, ?, ?, ?> getUniqueExtension(
			UnivariatePolynomial<T> minimalPolynomial) {
		throw new UnsupportedOperationException("Not implemented!");
	}

	@Override
	public T getRandomInteger() {
		return getRandomElement();
	}

	@Override
	public DiscreteValuationRing<T, T> ringOfIntegers() {
		return new LocalRingImplementation<>(this, field.toString());
	}

}
