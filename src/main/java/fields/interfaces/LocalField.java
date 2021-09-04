package fields.interfaces;

import fields.local.Value;

public interface LocalField<T extends Element<T>, S extends Element<S>> extends ValueField<T> {
	boolean isComplete();

	T inverse(T t, int accuracy);

	T negative(T t, int accuracy);

	boolean isInteger(T t);

	public interface Valuation<T extends Element<T>> {
		Value valuation(T t);
	}

	Value valuation(T t);

	default Valuation<T> valuation() {
		return new Valuation<>() {
			@Override
			public Value valuation(T t) {
				return LocalField.this.valuation(t);
			}
		};
	}

	T uniformizer();

	public Field<S> reduction();

	public S reduce(T t);
	
	public T upToUniformizerPower(T t);

	public T lift(S s);

	public T round(T t, int accuracy);

	public int getAccuracy();

	public LocalField<T, S> withAccuracy(int accuracy);

	public class OtherVersion<T extends Element<T>, U extends Element<U>, S extends Element<S>, LF extends LocalField<U, S>> {
		private LF field;
		private MathMap<T, U> retraction;
		private MathMap<U, T> embedding;

		public OtherVersion(LF field, MathMap<T, U> retraction, MathMap<U, T> embedding) {
			this.field = field;
			this.retraction = retraction;
			this.embedding = embedding;
		}

		public LF getField() {
			return field;
		}

		public MathMap<T, U> getRetraction() {
			return retraction;
		}

		public MathMap<U, T> getEmbedding() {
			return embedding;
		}
	}

	public OtherVersion<T, ?, S, ?> exact();

	public OtherVersion<T, ?, S, ?> complete(int accuracy);

	T getRandomInteger();

	LocalRing<T, S> ringOfIntegers();
	
}
