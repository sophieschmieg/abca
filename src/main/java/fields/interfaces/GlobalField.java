package fields.interfaces;

import java.util.ArrayList;
import java.util.List;

import fields.helper.FieldEmbedding;
import fields.vectors.Vector;

public interface GlobalField<T extends Element<T>, I extends Element<I>, S extends Element<S>> extends Field<T> {
	T getInteger(I t);

	DedekindRing<I, T, S> ringOfIntegers();

	public static class ExtensionOfGlobalField<T extends Element<T>, IT extends Element<IT>, ST extends Element<ST>, B extends Element<B>, I extends Element<I>, S extends Element<S>, E extends AlgebraicExtensionElement<B, E>, IE extends Element<IE>, R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>, DFE extends DiscreteValuationFieldExtension<B, S, E, DFE, R, RE, RFE>, FE extends GlobalFieldExtension<B, I, S, E, IE, R, RE, RFE, DFE, FE>> {
		private GlobalField<T, IT, ST> field;
		private FE extension;
		private MathMap<T, E> embeddingMap;
		private MathMap<E, Vector<T>> asVectorMap;

		public ExtensionOfGlobalField(GlobalField<T, IT, ST> field, FE extension, MathMap<T, E> embeddingMap,
				MathMap<E, Vector<T>> asVectorMap) {
			this.field = field;
			this.extension = extension;
			this.embeddingMap = embeddingMap;
			this.asVectorMap = asVectorMap;
		}

		public GlobalField<T, IT, ST> getField() {
			return field;
		}

		public FE getExtension() {
			return extension;
		}

		public MathMap<T, E> getEmbeddingMap() {
			return embeddingMap;
		}
		
		public E getEmbedding(T t) {
			return embeddingMap.evaluate(t);
		}

		public MathMap<E, Vector<T>> getAsVectorMap() {
			return asVectorMap;
		}
		
		public Vector<T> asVector(E t) {
			return asVectorMap.evaluate(t);
		}

		public Extension<T, B, E, FE> asExtension() {
			return new Extension<>(extension, field, embeddingMap, asVectorMap);
		}

		public ExtensionOfGlobalField<T, IT, ST, B, I, S, E, IE, R, RE, RFE, DFE, FE> extendFurther(
				UnivariatePolynomial<E> minimalPolynomial) {
			return extendFurther(extension.getEmbeddedExtension(minimalPolynomial));
		}
		public ExtensionOfGlobalField<T, IT, ST, B, I, S, E, IE, R, RE, RFE, DFE, FE> extendFurther(
				FieldEmbedding<B, E, FE> embedding) {
			return new ExtensionOfGlobalField<>(field, embedding.getField(), new MathMap<>() {
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

	ExtensionOfGlobalField<T, I, S, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?> getGlobalFieldExtension(
			UnivariatePolynomial<T> minimalPolynomial);
}
