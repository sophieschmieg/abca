package fields.interfaces;

import fields.local.Value;
import fields.vectors.Vector;

public interface DiscreteValuationField<T extends Element<T>, S extends Element<S>> extends ValueField<T> {
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
				return DiscreteValuationField.this.valuation(t);
			}
		};
	}

	T uniformizer();

	public Field<S> residueField();

	public S reduceInteger(T t);
	
	public T upToUniformizerPower(T t);

	public T liftToInteger(S s);

	public T round(T t, int accuracy);

	public int getAccuracy();

	public DiscreteValuationField<T, S> withAccuracy(int accuracy);

	public class OtherVersion<T extends Element<T>, U extends Element<U>, S extends Element<S>, LF extends DiscreteValuationField<U, S>> {
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
	
	public class DiscreteValuationFieldExtension<T extends Element<T>, S extends Element<S>, B extends Element<B>, E extends AlgebraicExtensionElement<B, E>, RE extends Element<RE>, FE extends FieldExtension<B, E, FE> & DiscreteValuationField<E, RE>> {
		private DiscreteValuationField<T, S> field;
		private FE extension;
		private MathMap<T, E> embeddingMap;
		private MathMap<E, Vector<T>> asVectorMap;
		
		public DiscreteValuationFieldExtension(DiscreteValuationField<T, S> field, FE extension,
				MathMap<T, E> embeddingMap, MathMap<E, Vector<T>> asVectorMap) {
			this.field = field;
			this.extension = extension;
			this.embeddingMap = embeddingMap;
			this.asVectorMap = asVectorMap;
		}

		public DiscreteValuationField<T, S> getField() {
			return field;
		}

		public FE getExtension() {
			return extension;
		}

		public MathMap<T, E> getEmbeddingMap() {
			return embeddingMap;
		}

		public MathMap<E, Vector<T>> getAsVectorMap() {
			return asVectorMap;
		}
		
	}
	
	DiscreteValuationFieldExtension<T, S, ?, ?, ?, ?> getUniqueExtension(UnivariatePolynomial<T> minimalPolynomial);
	
	public OtherVersion<T, ?, S, ?> exact();

	public OtherVersion<T, ?, S, ?> complete(int accuracy);

	T getRandomInteger();

	DiscreteValuationRing<T, S> ringOfIntegers();
	
}
