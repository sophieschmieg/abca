package fields.helper;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.FieldOfFractions.Fraction;
import fields.helper.TranscendentalFieldExtension.TExt;
import fields.interfaces.Algebra;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.VectorSpace;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;
import fields.vectors.MatrixAlgebra;
import fields.vectors.Vector;

public class TranscendentalFieldExtension<T extends Element<T>> extends AbstractField<TExt<T>>
		implements Algebra<T, TExt<T>>, Field<TExt<T>>, VectorSpace<T, TExt<T>> {
	public static class TExt<T extends Element<T>> extends AbstractElement<TExt<T>> {
		private Fraction<Polynomial<T>> asFraction;

		private TExt(Fraction<Polynomial<T>> asFraction) {
			this.asFraction = asFraction;
		}

		@SuppressWarnings("unchecked")
		@Override
		public boolean equals(Object o) {
			return asFraction.equals(((TExt<T>) o).asFraction);
		}

		@Override
		public int compareTo(TExt<T> o) {
			return asFraction.compareTo(o.asFraction);
		}

		@Override
		public String toString() {
			return asFraction.toString();
		}
		
		public Polynomial<T> getNumerator() {
			return asFraction.getNumerator();
		}
		
		public Polynomial<T> getDenominator() {
			return asFraction.getDenominator();
		}
		
		public Polynomial<T> asInteger() {
			return asFraction.asInteger();
		}
	}

	private FieldOfFractions<Polynomial<T>> asFieldOfFractions;
	private Field<T> baseField;
	private PolynomialRing<T> polynomialRing;
	private TExt<T> zero;
	private TExt<T> one;

	public TranscendentalFieldExtension(Field<T> baseField, int dimensions) {
		this(baseField, dimensions, Monomial.GREVLEX);
	}

	public TranscendentalFieldExtension(Field<T> baseField, int dimensions, Comparator<Monomial> comparator) {
		this(baseField, AbstractPolynomialRing.getPolynomialRing(baseField, dimensions, comparator));
	}

	public TranscendentalFieldExtension(Field<T> baseField, PolynomialRing<T> polynomialRing) {
		this.baseField = baseField;
		this.polynomialRing = polynomialRing;
		this.asFieldOfFractions = new FieldOfFractions<>(polynomialRing);
		this.zero = new TExt<>(asFieldOfFractions.zero());
		this.one = new TExt<>(asFieldOfFractions.one());
	}

	@Override
	public Field<T> getRing() {
		return getField();
	}

	@Override
	public TExt<T> zero() {
		return zero;
	}

	@Override
	public TExt<T> add(TExt<T> s1, TExt<T> s2) {
		return new TExt<>(asFieldOfFractions.add(s1.asFraction, s2.asFraction));
	}

	@Override
	public TExt<T> negative(TExt<T> s) {
		return new TExt<>(asFieldOfFractions.negative(s.asFraction));
	}

	@Override
	public TExt<T> scalarMultiply(T t, TExt<T> s) {
		return new TExt<>(asFieldOfFractions.multiply(getEmbedding(t).asFraction, s.asFraction));
	}

	@Override
	public boolean isFree() {
		return true;
	}

	@Override
	public TExt<T> scalarMultiply(int n, TExt<T> s) {
		return new TExt<>(asFieldOfFractions.multiply(getInteger(n).asFraction, s.asFraction));
	}

	@Override
	public TExt<T> scalarMultiply(BigInteger n, TExt<T> s) {
		return new TExt<>(asFieldOfFractions.multiply(getInteger(n).asFraction, s.asFraction));
	}

	@Override
	public TExt<T> scalarMultiply(T t1, T t2, TExt<T> s) {
		return new TExt<>(
				asFieldOfFractions.multiply(getEmbedding(baseField.multiply(t1, t2)).asFraction, s.asFraction));
	}

	@Override
	public TExt<T> scalarMultiply(int n, T t, TExt<T> s) {
		return new TExt<>(asFieldOfFractions.multiply(getEmbedding(baseField.multiply(n, t)).asFraction, s.asFraction));
	}

	@Override
	public TExt<T> scalarMultiply(BigInteger n, T t, TExt<T> s) {
		return new TExt<>(asFieldOfFractions.multiply(getEmbedding(baseField.multiply(n, t)).asFraction, s.asFraction));
	}

	@Override
	public TExt<T> scalarMultiply(int n, T t1, T t2, TExt<T> s) {
		return new TExt<>(
				asFieldOfFractions.multiply(getEmbedding(baseField.multiply(n, t1, t2)).asFraction, s.asFraction));
	}

	@Override
	public TExt<T> scalarMultiply(BigInteger n, T t1, T t2, TExt<T> s) {
		return new TExt<>(
				asFieldOfFractions.multiply(getEmbedding(baseField.multiply(n, t1, t2)).asFraction, s.asFraction));
	}

	@Override
	public boolean isLinearIndependent(List<TExt<T>> s) {
		Polynomial<T> lcm = polynomialRing.one();
		for (TExt<T> element : s) {
			lcm = polynomialRing.lcm(lcm, element.asFraction.getDenominator());
		}
		List<Polynomial<T>> numerators = new ArrayList<>();
		for (TExt<T> element : s) {
			numerators.add(
					asFieldOfFractions.multiply(asFieldOfFractions.getEmbedding(lcm), element.asFraction).asInteger());
		}
		return polynomialRing.isLinearIndependent(numerators);
	}

	@Override
	public boolean isGeneratingModule(List<TExt<T>> s) {
		return false;
	}

	@Override
	public List<T> nonTrivialCombination(List<TExt<T>> s) {
		return nonTrivialCombinations(s).get(0);
	}

	@Override
	public List<List<T>> nonTrivialCombinations(List<TExt<T>> s) {
		Polynomial<T> lcm = polynomialRing.one();
		for (TExt<T> element : s) {
			lcm = polynomialRing.lcm(lcm, element.asFraction.getDenominator());
		}
		List<Polynomial<T>> numerators = new ArrayList<>();
		for (TExt<T> element : s) {
			numerators.add(
					asFieldOfFractions.multiply(asFieldOfFractions.getEmbedding(lcm), element.asFraction).asInteger());
		}
		return polynomialRing.nonTrivialCombinations(numerators);
	}

	@Override
	public List<TExt<T>> getModuleGenerators() {
		throw new InfinityException();
	}

	@Override
	public Vector<T> asVector(TExt<T> s) {
		throw new UnsupportedOperationException("basis only exists due to axiom of choice");
	}

	@Override
	public TExt<T> fromVector(Vector<T> vector) {
		throw new UnsupportedOperationException("basis only exists due to axiom of choice");
	}

	@Override
	public Exactness exactness() {
		return baseField.exactness();
	}

	@Override
	public TExt<T> getRandomElement() {
		return new TExt<>(asFieldOfFractions.getRandomElement());
	}

	@Override
	public boolean isFinite() {
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		throw new InfinityException();
	}

	@Override
	public Iterator<TExt<T>> iterator() {
		throw new InfinityException();
	}

	@Override
	public TExt<T> one() {
		return one;
	}

	@Override
	public BigInteger characteristic() {
		return baseField.characteristic();
	}

	@Override
	public TExt<T> multiply(TExt<T> t1, TExt<T> t2) {
		return new TExt<>(asFieldOfFractions.multiply(t1.asFraction, t2.asFraction));
	}

	@Override
	public TExt<T> inverse(TExt<T> t) {
		return new TExt<>(asFieldOfFractions.inverse(t.asFraction));
	}

	@Override
	public TExt<T> getEmbedding(T t) {
		return new TExt<>(asFieldOfFractions.getEmbedding(polynomialRing.getEmbedding(t)));
	}
	
	public TExt<T> getEmbedding(Polynomial<T> t) {
		return new TExt<>(asFieldOfFractions.getEmbedding(t));
	}

	@Override
	public boolean isGeneratingAlgebra(List<TExt<T>> s) {
		return false;
	}

	@Override
	public List<TExt<T>> getAlgebraGenerators() {
		throw new InfinityException();
	}

	@Override
	public Field<T> getField() {
		return baseField;
	}

	@Override
	public List<TExt<T>> getBasis() {
		throw new InfinityException();
	}

	@Override
	public int dimension() {
		throw new InfinityException();
	}

	public int transcendenceDegree() {
		return polynomialRing.numberOfVariables();
	}
	
	public PolynomialRing<T> polynomialRing() {
		return polynomialRing;
	}

	@Override
	public MatrixAlgebra<T> matrixAlgebra() {
		throw new InfinityException();
	}
}
