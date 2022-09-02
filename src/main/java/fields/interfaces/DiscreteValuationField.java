package fields.interfaces;

import java.util.ArrayList;
import java.util.List;

import fields.helper.FieldEmbedding;
import fields.local.Value;
import fields.vectors.Vector;

public interface DiscreteValuationField<T extends Element<T>, S extends Element<S>>
		extends ValueField<T>, GlobalField<T, T, S> {
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

	public class ExtensionOfDiscreteValuationField<T extends Element<T>, SR extends Element<SR>, B extends Element<B>, S extends Element<S>, E extends AlgebraicExtensionElement<B, E>, FE extends DiscreteValuationFieldExtension<B, S, E, FE, R, RE, RFE>, R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>> {
		private DiscreteValuationField<T, SR> field;
		private FE extension;
		private MathMap<T, E> embeddingMap;
		private MathMap<E, Vector<T>> asVectorMap;

		public ExtensionOfDiscreteValuationField(DiscreteValuationField<T, SR> field, FE extension,
				MathMap<T, E> embeddingMap, MathMap<E, Vector<T>> asVectorMap) {
			this.field = field;
			this.extension = extension;
			this.embeddingMap = embeddingMap;
			this.asVectorMap = asVectorMap;
		}

		public DiscreteValuationField<T, SR> getField() {
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

		public ExtensionOfGlobalField<T, T, SR, B, B, S, E, E, R, RE, RFE, FE, FE> asGlobalFieldExtension() {
			return new ExtensionOfGlobalField<T, T, SR, B, B, S, E, E, R, RE, RFE, FE, FE>(field, extension,
					embeddingMap, asVectorMap);
		}

		public Extension<T, B, E, FE> asExtension() {
			return new Extension<>(extension, field, embeddingMap, asVectorMap);
		}

		public ExtensionOfDiscreteValuationField<T, SR, B, S, E, FE, R, RE, RFE> extendFurther(
				FieldEmbedding<B, E, FE> embedding) {
			return new ExtensionOfDiscreteValuationField<>(field, embedding.getField(), new MathMap<>() {
				@Override
				public E evaluate(T t) {
					return embedding.getEmbedding(embeddingMap.evaluate(t));
				}
			}, new MathMap<>() {
				@Override
				public Vector<T> evaluate(E t) {
					Vector<E> asVector = embedding.asVector(t);
					List<T> result = new ArrayList<>();
					for (E element : asVector.asList()) {
						result.addAll(asVectorMap.evaluate(element).asList());
					}
					return new Vector<>(result);
				}
			});
		}
	}

	ExtensionOfDiscreteValuationField<T, S, ?, ?, ?, ?, ?, ?, ?> getUniqueExtension(
			UnivariatePolynomial<T> minimalPolynomial);

	@Override
	default ExtensionOfGlobalField<T, T, S, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?> getGlobalFieldExtension(
			UnivariatePolynomial<T> minimalPolynomial) {
		return getUniqueExtension(minimalPolynomial).asGlobalFieldExtension();
	}

	public OtherVersion<T, ?, S, ?> exact();

	public OtherVersion<T, ?, S, ?> complete(int accuracy);

	T getRandomInteger();

	DiscreteValuationRing<T, S> ringOfIntegers();

	@Override
	default T getInteger(T t) {
		return t;
	}

}
