package fields.interfaces;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import fields.exceptions.InfinityException;
import fields.helper.FieldEmbedding;
import fields.integers.Integers.IntE;
import fields.integers.Rationals.Fraction;
import fields.vectors.Vector;
import util.ConCatMap;

public interface Field<T extends Element<T>> extends MathSet<T>, Ring<T>, DedekindRing<T, T, T>, LocalRing<T, T, T> {
	public T zero();

	public T one();

	public BigInteger characteristic();

	public T add(T t1, T t2);

	public T negative(T t);

	public T multiply(T t1, T t2);

	public T inverse(T t);

	public Group<T> getAdditiveGroup();

	public Group<T> getMultiplicativeGroup();

	public T add(T t1, T t2, T t3);

	public T add(T t1, T t2, T t3, T t4);

	public T subtract(T minuend, T subtrahend);

	public T getInteger(int n);

	public T getInteger(BigInteger n);

	public T getInteger(IntE n);

	public T getFraction(Fraction t);

	public T multiply(int n, T t);

	public T multiply(BigInteger n, T t);

	public T multiply(IntE n, T t);

	public T multiply(T t1, T t2, T t3);

	public T multiply(int n, T t1, T t2);

	public T multiply(BigInteger n, T t1, T t2);

	public T multiply(IntE n, T t1, T t2);

	public T multiply(int n, T t1, T t2, T t3);

	public T multiply(BigInteger n, T t1, T t2, T t3);

	public T multiply(IntE n, T t1, T t2, T t3);

	public T divide(T dividend, T divisor);

	public T power(T t, int n);

	public T power(T t, BigInteger n);

	public Iterable<T> getNonZeroElements() throws InfinityException;

	public boolean isIrreducible(UnivariatePolynomial<T> t);

	public FactorizationResult<Polynomial<T>, T> factorization(UnivariatePolynomial<T> t);

	public boolean hasRoots(Polynomial<T> t);

	public Map<T, Integer> roots(Polynomial<T> t);

	public boolean hasRoot(T t, int n);

	public Map<T, Integer> roots(T t, int n);

	public boolean hasSqrt(T t);

	public Map<T, Integer> sqrt(T t);

	public T primitiveRoot();

	public T characteristicRoot(T t);

	public static class Extension<T extends Element<T>, Base extends Element<Base>, Ext extends AlgebraicExtensionElement<Base, Ext>, ExtField extends FieldExtension<Base, Ext, ExtField>> {
		private ExtField extension;
		private Field<T> domain;
		private MathMap<T, Ext> embeddingMap;
		private MathMap<Ext, Vector<T>> asVectorMap;
		private MathMap<Ext, T> retractionMap;

		public Extension(ExtField extension, Field<T> domain, MathMap<T, Ext> embeddingMap,
				MathMap<Ext, Vector<T>> asVectorMap) {
			this.extension = extension;
			this.domain = domain;
			this.embeddingMap = embeddingMap;
			this.asVectorMap = asVectorMap;
			this.retractionMap = new ConCatMap<Ext, Vector<T>, T>(asVectorMap, new MathMap<>() {
				@Override
				public T evaluate(Vector<T> t) {
					return t.get(1);
				}
			});
		}

		public ExtField extension() {
			return extension;
		}

		public MathMap<T, Ext> embeddingMap() {
			return embeddingMap;
		}

		public MathMap<Ext, Vector<T>> asVectorMap() {
			return asVectorMap;
		}

		public MathMap<Ext, T> retractionMap() {
			return retractionMap;
		}

		@Override
		public String toString() {
			return extension.toString() + "/" + domain.toString();
		}

		public Field<T> domain() {
			return domain;
		}

		public Extension<T, Base, Ext, ExtField> extendFurther(FieldEmbedding<Base, Ext, ExtField> embedding) {
			return new Extension<>(embedding.getField(), domain, new MathMap<>() {
				@Override
				public Ext evaluate(T t) {
					return embedding.getEmbedding(embeddingMap.evaluate(t));
				}
			}, new MathMap<>() {
				@Override
				public Vector<T> evaluate(Ext t) {
					Vector<Ext> asVector = embedding.asVector(t);
					List<T> result = new ArrayList<>();
					for (Ext element : asVector.asList()) {
						result.addAll(asVectorMap.evaluate(element).asList());
					}
					return new Vector<>(result);
				}
			});
		}

	}

	public Extension<T, ?, ?, ?> getExtension(UnivariatePolynomial<T> minimalPolynomial);
}
