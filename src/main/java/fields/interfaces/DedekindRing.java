package fields.interfaces;

import fields.local.Value;

public interface DedekindRing<T extends Element<T>, U extends Element<U>, S extends Element<S>> extends Ring<T> {
	Value valuation(T t, Ideal<T> maximalIdeal);

	Field<U> fieldOfFractions();
	
	MathMap<T, U> embedding();
	boolean isInteger(U t);
	T asInteger(U t);

	Field<S> reduction(Ideal<T> maximalIdeal);

	S reduce(T t, Ideal<T> maximalIdeal);

	T lift(S s, Ideal<T> maximalIdeal);
	
	LocalRing<U, S> localize(Ideal<T> maximalIdeal);
}
