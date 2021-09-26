package fields.interfaces;

import fields.local.Value;

public interface DedekindRing<T extends Element<T>, U extends Element<U>, S extends Element<S>> extends Ring<T> {
	Value valuation(T t, Ideal<T> maximalIdeal);

	Field<U> fieldOfFractions();

	MathMap<T, U> embedding();

	boolean isInteger(U t);
	
	@Override
	default boolean isIntegral() {
		return true;
	}

	@Override
	default boolean isReduced() {
		return true;
	}

	@Override
	default boolean isIrreducible() {
		return true;
	}
	
	@Override
	default boolean isDedekindDomain() {
		return true;
	}

	T asInteger(U t);

	Field<S> reduction(Ideal<T> maximalIdeal);

	S reduce(T t, Ideal<T> maximalIdeal);

	T lift(S s, Ideal<T> maximalIdeal);

	LocalRing<U, S> localize(Ideal<T> maximalIdeal);
}
