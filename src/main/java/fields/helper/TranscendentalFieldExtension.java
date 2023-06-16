package fields.helper;

import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;

import fields.exceptions.InfinityException;
import fields.helper.FieldOfFractions.Fraction;
import fields.helper.TranscendentalFieldExtension.TExt;
import fields.interfaces.Algebra;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Ideal;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.interfaces.VectorSpace;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;
import util.PeekableReader;
import varieties.SimpleFunctionField;
import varieties.SimpleFunctionField.SimpleRationalFunction;

public class TranscendentalFieldExtension<T extends Element<T>> extends AbstractField<TExt<T>>
		implements Algebra<T, TExt<T>>, Field<TExt<T>>, VectorSpace<T, TExt<T>> {
	public static class TExt<T extends Element<T>> extends AbstractElement<TExt<T>> {
		private Fraction<Polynomial<T>> asFraction;

		private TExt(TranscendentalFieldExtension<T> extension, Fraction<Polynomial<T>> asFraction) {
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

	public TranscendentalFieldExtension(Field<T> baseField, String[] variables) {
		this(baseField, variables, Monomial.GREVLEX);
	}

	public TranscendentalFieldExtension(Field<T> baseField, int dimensions) {
		this(baseField, dimensions, Monomial.GREVLEX);
	}

	public TranscendentalFieldExtension(Field<T> baseField, String[] variables, Comparator<Monomial> comparator) {
		this(baseField, AbstractPolynomialRing.getPolynomialRing(baseField, comparator, variables));
	}

	public TranscendentalFieldExtension(Field<T> baseField, int dimensions, Comparator<Monomial> comparator) {
		this(baseField, AbstractPolynomialRing.getPolynomialRing(baseField, dimensions, comparator));
	}

	public TranscendentalFieldExtension(Field<T> baseField, PolynomialRing<T> polynomialRing) {
		this.baseField = baseField;
		this.polynomialRing = polynomialRing;
		this.asFieldOfFractions = new FieldOfFractions<>(polynomialRing);
		this.zero = getElement(asFieldOfFractions.zero());
		this.one = getElement(asFieldOfFractions.one());
	}

	@Override
	public String toString() {
		return baseField + "(" + String.join(",", polynomialRing.getVariableNames()) + ")";
	}

	@Override
	public TExt<T> parse(PeekableReader reader) throws IOException {
		return getElement(asFieldOfFractions.parse(reader));
	}

	@Override
	public Field<T> getRing() {
		return getField();
	}

	@Override
	public TExt<T> zero() {
		return zero;
	}

	private TExt<T> getElement(Fraction<Polynomial<T>> asFraction) {
		Polynomial<T> numerator = asFraction.getNumeratorDirect();
		Polynomial<T> denominator = asFraction.getDenominatorDirect();
		T lc = denominator.leadingCoefficient();
		numerator = polynomialRing.divideScalar(numerator, lc);
		denominator = polynomialRing.normalize(denominator);
		return new TExt<>(this, asFieldOfFractions.getFraction(numerator, denominator));
	}

	@Override
	public TExt<T> add(TExt<T> s1, TExt<T> s2) {
		return getElement(asFieldOfFractions.add(s1.asFraction, s2.asFraction));
	}

	@Override
	public TExt<T> negative(TExt<T> s) {
		return getElement(asFieldOfFractions.negative(s.asFraction));
	}

	@Override
	public TExt<T> scalarMultiply(T t, TExt<T> s) {
		return getElement(asFieldOfFractions.multiply(getEmbedding(t).asFraction, s.asFraction));
	}

	@Override
	public boolean isFree() {
		return true;
	}

	@Override
	public boolean isTorsionFree() {
		return true;
	}

	@Override
	public Ideal<T> annihilator() {
		return baseField.getZeroIdeal();
	}

	@Override
	public TExt<T> scalarMultiply(int n, TExt<T> s) {
		return getElement(asFieldOfFractions.multiply(getInteger(n).asFraction, s.asFraction));
	}

	@Override
	public TExt<T> scalarMultiply(BigInteger n, TExt<T> s) {
		return getElement(asFieldOfFractions.multiply(getInteger(n).asFraction, s.asFraction));
	}

	@Override
	public TExt<T> scalarMultiply(T t1, T t2, TExt<T> s) {
		return getElement(
				asFieldOfFractions.multiply(getEmbedding(baseField.multiply(t1, t2)).asFraction, s.asFraction));
	}

	@Override
	public TExt<T> scalarMultiply(int n, T t, TExt<T> s) {
		return getElement(asFieldOfFractions.multiply(getEmbedding(baseField.multiply(n, t)).asFraction, s.asFraction));
	}

	@Override
	public TExt<T> scalarMultiply(BigInteger n, T t, TExt<T> s) {
		return getElement(asFieldOfFractions.multiply(getEmbedding(baseField.multiply(n, t)).asFraction, s.asFraction));
	}

	@Override
	public TExt<T> scalarMultiply(int n, T t1, T t2, TExt<T> s) {
		return getElement(
				asFieldOfFractions.multiply(getEmbedding(baseField.multiply(n, t1, t2)).asFraction, s.asFraction));
	}

	@Override
	public TExt<T> scalarMultiply(BigInteger n, T t1, T t2, TExt<T> s) {
		return getElement(
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

//	@Override
//	public List<T> nonTrivialCombination(List<TExt<T>> s) {
//		return nonTrivialCombinations(s).get(0);
//	}

	@Override
	public List<Vector<T>> nonTrivialCombinations(List<TExt<T>> s) {
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
		return getElement(asFieldOfFractions.getRandomElement());
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
		return getElement(asFieldOfFractions.multiply(t1.asFraction, t2.asFraction));
	}

	@Override
	public TExt<T> inverse(TExt<T> t) {
		return getElement(asFieldOfFractions.inverse(t.asFraction));
	}

	@Override
	public TExt<T> getEmbedding(T t) {
		return getElement(asFieldOfFractions.getEmbedding(polynomialRing.getEmbedding(t)));
	}

	public MathMap<T, TExt<T>> getEmbeddingMap() {
		return new MathMap<>() {
			@Override
			public TExt<T> evaluate(T t) {
				return getEmbedding(t);
			}
		};
	}

	public TExt<T> getEmbedding(Polynomial<T> t) {
		return getElement(asFieldOfFractions.getEmbedding(t));
	}

	public TExt<T> getElement(Polynomial<T> numerator, Polynomial<T> denomiator) {
		return getElement(asFieldOfFractions.getFraction(numerator, denomiator));
	}

	public <S extends Element<S>> TExt<T> getEmbedding(TExt<S> t, MathMap<S, T> map) {
		return getElement(polynomialRing.getEmbedding(t.getNumerator(), map),
				polynomialRing.getEmbedding(t.getDenominator(), map));
	}

	public TExt<T> getVar(int i) {
		return getEmbedding(polynomialRing.getVar(i));
	}

	public TExt<T> getVarPower(int i, int power) {
		if (power < 0) {
			return inverse(getVarPower(i, -power));
		}
		return getEmbedding(polynomialRing.getVarPower(i, power));
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
	public TExt<T> getUnitVector(int index) {
		throw new InfinityException();
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

	@Override
	public Extension<TExt<T>, TExt<T>, SimpleRationalFunction<T>, SimpleFunctionField<T>> getExtension(
			UnivariatePolynomial<TExt<T>> minimalPolynomial) {
		SimpleFunctionField<T> extension = new SimpleFunctionField<>(minimalPolynomial, this, "x");
		return new Extension<>(extension, this, extension.getEmbeddingMap(), extension.asVectorMap());
	}

	@Override
	public FactorizationResult<Polynomial<TExt<T>>, TExt<T>> factorization(UnivariatePolynomial<TExt<T>> t) {
		UnivariatePolynomialRing<Fraction<Polynomial<T>>> univariate = asFieldOfFractions.getUnivariatePolynomialRing();
		UnivariatePolynomial<Fraction<Polynomial<T>>> polynomial = univariate.getEmbedding(t, new MathMap<>() {
			@Override
			public Fraction<Polynomial<T>> evaluate(TExt<T> t) {
				return t.asFraction;
			}
		});
		FactorizationResult<Polynomial<Fraction<Polynomial<T>>>, Fraction<Polynomial<T>>> factorization = asFieldOfFractions
				.factorization(polynomial);
		SortedMap<Polynomial<TExt<T>>, Integer> result = new TreeMap<>();
		for (Polynomial<Fraction<Polynomial<T>>> factor : factorization.primeFactors()) {
			result.put(getUnivariatePolynomialRing().getEmbedding(factor, new MathMap<>() {
				@Override
				public TExt<T> evaluate(Fraction<Polynomial<T>> t) {
					return getElement(t);
				}
			}), factorization.multiplicity(factor));
		}
		return new FactorizationResult<>(getElement(factorization.getUnit()), result);
	}

	public TExt<T> characteristicRoot(TExt<T> t, int power) {
		return divide(getEmbedding(polynomialRing.characteristicRoot(t.getNumerator(), power)),
				getEmbedding(polynomialRing.characteristicRoot(t.getDenominator(), power)));
	}

	@Override
	public boolean hasCharacteristicRoot(TExt<T> t, int power) {
		return polynomialRing.hasCharacteristicRoot(t.getNumerator(), power)
				&& polynomialRing.hasCharacteristicRoot(t.getDenominator(), power);
	}

	public TExt<T> substitute(Polynomial<T> t, List<TExt<T>> substitutes) {
		TExt<T> result = zero();
		for (Monomial m : t.monomials()) {
			TExt<T> value = getEmbedding(t.coefficient(m));
			for (int i = 0; i < m.exponents().length; i++) {
				if (m.exponents()[i] != 0) {
					value = multiply(power(substitutes.get(i), m.exponents()[i]), value);
				}
			}
			result = add(value, result);
		}
		return result;
	}

	public TExt<T> substitute(TExt<T> t, List<TExt<T>> substitutes) {
		return divide(substitute(t.getNumerator(), substitutes), substitute(t.getDenominator(), substitutes));
	}

	private Matrix<T> asMatrix(List<TExt<T>> m) {
		List<Vector<T>> result = new ArrayList<>();
		for (TExt<T> s : m) {
			result.add(asVector(s));
		}
		return Matrix.fromColumns(result);
	}
	
	@Override
	public boolean isSubModuleMemberModule(List<TExt<T>> m, TExt<T> b) {
		Matrix<T> matrix = asMatrix(m);
		return isSubModuleMemberModule(matrix.getModule(getRing()), matrix, asVector(b));
	}
	
	@Override
	public boolean isSubModuleMemberModule(MatrixModule<T> module, Matrix<T> m, Vector<T> b) {
		return getRing().isSubModuleMember(module, m, b);
	}
	
	@Override
	public TExt<T> asSubModuleMemberModule(List<TExt<T>> m, TExt<T> b) {
		Matrix<T> matrix = asMatrix(m);
		return fromVector(asSubModuleMemberModule(matrix.getModule(getRing()), matrix, asVector(b)));
	}
	
	@Override
	public Vector<T> asSubModuleMemberModule(MatrixModule<T> module, Matrix<T> m, Vector<T> b) {
		return getRing().asSubModuleMember(module, m, b);
	}
	
	@Override
	public List<Vector<T>> syzygyProblemModule(List<TExt<T>> m) {
		Matrix<T> matrix = asMatrix(m);
		return syzygyProblemModule(matrix.getModule(getRing()), matrix);
	}
	
	@Override
	public List<Vector<T>> syzygyProblemModule(MatrixModule<T> module, Matrix<T> m) {
		return getRing().syzygyProblem(module, m);
	}
	
	@Override
	public List<TExt<T>> simplifySubModuleGeneratorsModule(List<TExt<T>> m) {
		Matrix<T> matrix = asMatrix(m);
		List<Vector<T>> simplified = simplifySubModuleGeneratorsModule(matrix.getModule(getRing()), matrix);
		List<TExt<T>> result = new ArrayList<>();
		for (Vector<T> vector : simplified) {
			result.add(fromVector(vector));
		}
		return result;
	}
	
	@Override
	public List<Vector<T>> simplifySubModuleGeneratorsModule(MatrixModule<T> module, Matrix<T> m) {
		return getRing().simplifySubModuleGenerators(module, m);
	}
}
