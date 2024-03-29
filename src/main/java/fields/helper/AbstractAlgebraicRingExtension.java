package fields.helper;

import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import fields.exceptions.InfinityException;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.AlgebraicRingExtension;
import fields.interfaces.Element;
import fields.interfaces.Ideal;
import fields.interfaces.MathMap;
import fields.interfaces.Module;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.interfaces.UnivariatePolynomialRing.ExtendedResultantResult;
import fields.local.Value;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.CoordinateRing;
import fields.polynomials.CoordinateRing.CoordinateIdeal;
import fields.polynomials.CoordinateRing.CoordinateRingElement;
import fields.polynomials.GenericUnivariatePolynomialRing;
import fields.polynomials.Monomial;
import fields.polynomials.PolynomialIdeal;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;
import util.PeekableReader;

public abstract class AbstractAlgebraicRingExtension<T extends Element<T>, S extends AlgebraicExtensionElement<T, S>, Ext extends AlgebraicRingExtension<T, S, Ext>>
		extends AbstractAlgebra<T, S> implements AlgebraicRingExtension<T, S, Ext> {
	private Ring<T> baseRing;
	private S zero;
	private S one;
	private S alpha;
	private UnivariatePolynomialRing<T> polynomials;
	private UnivariatePolynomial<T> minimalPolynomial;
	private int degree;
	private Map<Integer, S> powers;
	private FreeModule<T> asFreeModule;
	private MatrixAlgebra<T> algebra;
	private PolynomialRing<T> genericPolynomialRing;
	private Vector<Polynomial<T>> genericVector;
	private Matrix<Polynomial<T>> genericMatrix;
	private FreeModule<Polynomial<T>> genericFreeModule;
	private MathMap<T, S> embedding;
	private FactorizationResult<Polynomial<T>, T> factorizedMinimalPolynomial;
	private List<Ext> asIrreducibleProduct;
	private List<S> irreducibleProductInverses;
	private String variableName;

	public AbstractAlgebraicRingExtension(UnivariatePolynomial<T> minimalPolynomial, Ring<T> baseRing,
			String variableName) {
		this.baseRing = baseRing;
		this.polynomials = baseRing.getUnivariatePolynomialRing().withVariableName(variableName);
		this.variableName = variableName;
		this.minimalPolynomial = minimalPolynomial;
		this.degree = this.minimalPolynomial.degree();
		this.zero = fromSmallDegreePolynomial(polynomials.zero());
		this.one = fromSmallDegreePolynomial(polynomials.one());
		this.powers = new TreeMap<>();
		for (int i = 0; i < 2 * degree; i++) {
			this.powers.put(i, fromSmallDegreePolynomial(polynomials.toUnivariate(
					polynomials.reduce(polynomials.getVarPower(i), Collections.singletonList(minimalPolynomial)))));
		}
		this.alpha = fromPolynomial(polynomials.getVar());
		this.asFreeModule = new FreeModule<>(this.baseRing, degree);
		this.embedding = new MathMap<>() {
			@Override
			public S evaluate(T t) {
				return getEmbedding(t);
			}
		};
	}

	public AbstractAlgebraicRingExtension(Ring<T> baseRing) {
		this(baseRing.getUnivariatePolynomialRing().getVar(), baseRing, "X");
	}

	@Override
	public Exactness exactness() {
		return baseRing.exactness();
	}

	@Override
	public String getVariableName() {
		return variableName;
	}

	public PolynomialRing<T> genericPolynomialRing() {
		if (genericPolynomialRing == null) {
			this.genericPolynomialRing = AbstractPolynomialRing.getPolynomialRing(baseRing, degree(), Monomial.GREVLEX);
		}
		return genericPolynomialRing;
	}

	public FreeModule<Polynomial<T>> genericFreeModule() {
		if (genericFreeModule == null) {
			genericFreeModule = new FreeModule<>(genericPolynomialRing(), degree());
		}
		return genericFreeModule;
	}

	public Vector<Polynomial<T>> genericVector() {
		if (genericVector == null) {
			List<Polynomial<T>> asList = new ArrayList<>();
			for (int i = 0; i < degree(); i++) {
				asList.add(genericPolynomialRing.getVar(i + 1));
			}
			genericVector = new Vector<>(asList);
		}
		return genericVector;
	}

	public Matrix<Polynomial<T>> genericMatrix() {
		if (genericMatrix == null) {
			Vector<Polynomial<T>> genericVector = genericVector();
			UnivariatePolynomial<Polynomial<T>> genericMinimalPolynomial = genericPolynomialRing()
					.getUnivariatePolynomialRing().getEmbedding(minimalPolynomial(), new MathMap<>() {

						@Override
						public Polynomial<T> evaluate(T t) {
							return genericPolynomialRing().getEmbedding(t);
						}
					});
			genericMatrix = asMatrix(genericVector, genericMinimalPolynomial, genericPolynomialRing().zero(),
					genericFreeModule());
		}
		return genericMatrix;
	}

	@Override
	public final UnivariatePolynomial<T> minimalPolynomial() {
		return minimalPolynomial;
	}

	@Override
	public String toString() {
		if (degree() == 1) {
			return baseRing.toString();
		}
		return baseRing.getUnivariatePolynomialRing().toString() + "/(" + minimalPolynomial.toString() + ")";
	}

	@Override
	public S parse(PeekableReader reader) throws IOException {
		return fromPolynomial(polynomials.toUnivariate(polynomials.parse(reader)));
	}

	@Override
	public final FreeModule<T> asFreeModule() {
		return asFreeModule;
	}

	@Override
	public final MatrixAlgebra<T> matrixAlgebra() {
		if (algebra == null) {
			this.algebra = this.asFreeModule.matrixAlgebra();
		}
		return algebra;
	}

	@Override
	public final Ring<T> getRing() {
		return baseRing;
	}

	@Override
	public final boolean isFinite() {
		return baseRing.isFinite();
	}

	@Override
	public final BigInteger getNumberOfElements() throws InfinityException {
		return baseRing.getNumberOfElements().pow(degree);
	}

	@Override
	public final BigInteger characteristic() {
		return baseRing.characteristic();
	}

	@Override
	public S fromPolynomial(UnivariatePolynomial<T> polynomial) {
		if (polynomial.degree() < degree()) {
			return fromSmallDegreePolynomial(polynomial);
		}
		if (polynomial.degree() >= 2 * degree()) {
			return fromSmallDegreePolynomial(polynomials
					.toUnivariate(polynomials.reduce(polynomial, Collections.singletonList(minimalPolynomial))));
		}
		S result = fromSmallDegreePolynomial(polynomials.round(polynomial, degree));
		for (int i = degree; i <= polynomial.degree(); i++) {
			result = add(result, scalarMultiply(polynomial.univariateCoefficient(i), powers.get(i)));
		}
		return result;
	}

	protected abstract S fromSmallDegreePolynomial(UnivariatePolynomial<T> polynomial);

	@Override
	public final UnivariatePolynomial<T> asPolynomial(S s) {
		return s.asPolynomial();
	}

	@Override
	public final boolean isFree() {
		return true;
	}

	@Override
	public Ideal<T> annihilator() {
		return baseRing.getZeroIdeal();
	}

	@Override
	public final boolean isLinearIndependent(List<S> s) {
		List<Vector<T>> asVectors = new ArrayList<>();
		for (S element : s) {
			asVectors.add(asVector(element));
		}
		return asFreeModule.isLinearIndependent(asVectors);
	}

	@Override
	public final List<Vector<T>> nonTrivialCombinations(List<S> s) {
		List<Vector<T>> asVectors = new ArrayList<>();
		for (S t : s) {
			asVectors.add(asVector(t));
		}
		return asFreeModule.nonTrivialCombinations(asVectors);
	}

	@Override
	public List<Vector<T>> getSyzygies() {
		return Collections.emptyList();
	}

	@Override
	public final boolean isGeneratingModule(List<S> s) {
		List<Vector<T>> asVectors = new ArrayList<>();
		for (S element : s) {
			asVectors.add(asVector(element));
		}
		return asFreeModule.isGeneratingModule(asVectors);
	}

	@Override
	public final List<S> getModuleGenerators() {
		List<S> generators = new ArrayList<>();
		S s = one();
		for (int i = 0; i < degree; i++) {
			generators.add(s);
			s = multiply(s, alpha());
		}
		return generators;
	}

	@Override
	public final boolean isGeneratingAlgebra(List<S> s) {
		List<S> asModuleGenerators = new ArrayList<>();
		for (S element : s) {
			S generated = one();
			for (int i = 0; i < degree; i++) {
				asModuleGenerators.add(generated);
				generated = multiply(generated, element);
			}
		}
		return isGeneratingModule(asModuleGenerators);
	}

	@Override
	public final List<S> getAlgebraGenerators() {
		return Collections.singletonList(alpha());
	}

	@Override
	public final int degree() {
		return degree;
	}

	@Override
	public final Matrix<T> asMatrix(S t) {
		return asMatrix(asVector(t), minimalPolynomial(), baseRing.zero(), asFreeModule);
	}

	private static <Q extends Element<Q>> Matrix<Q> asMatrix(Vector<Q> t, UnivariatePolynomial<Q> minimalPolynomial,
			Q zero, Module<Q, Vector<Q>> asVectorSpace) {
		List<List<Q>> m = new ArrayList<>();
		Vector<Q> asVector = t;
		List<Q> minPolyList = new ArrayList<>();
		int degree = minimalPolynomial.degree();
		for (int i = 0; i < degree; i++) {
			m.add(new ArrayList<>());
			minPolyList.add(minimalPolynomial.univariateCoefficient(i));
		}
		Vector<Q> minPoly = new Vector<>(minPolyList);
		for (int j = 0; j < degree; j++) {
			List<Q> nextVector = new ArrayList<>();
			nextVector.add(zero);
			for (int i = 0; i < degree; i++) {
				m.get(i).add(asVector.get(i + 1));
				nextVector.add(asVector.get(i + 1));
			}
			Q high = nextVector.get(degree);
			asVector = asVectorSpace.subtract(new Vector<>(nextVector.subList(0, degree)),
					asVectorSpace.scalarMultiply(high, minPoly));
		}
		return new Matrix<Q>(m);
	}

	@Override
	public final MathMap<T, S> getEmbeddingMap() {
		return embedding;
	}

	@Override
	public final T asBaseFieldElement(S s) {
		if (s.asPolynomial().degree() > 0) {
			throw new ArithmeticException("not degree one!");
		}
		return s.asPolynomial().univariateCoefficient(0);
	}

	@Override
	public final MathMap<S, T> asBaseFieldElementMap() {
		return new MathMap<>() {

			@Override
			public T evaluate(S t) {
				return asBaseFieldElement(t);
			}
		};
	}

	@Override
	public MathMap<S, Vector<T>> asVectorMap() {
		return new MathMap<>() {
			@Override
			public Vector<T> evaluate(S t) {
				return asVector(t);
			}
		};
	}

	@Override
	public final UnivariatePolynomial<T> minimalPolynomial(S s) {
		return minimalPolynomial(s, degree, this, baseRing, new MathMap<>() {
			@Override
			public Vector<T> evaluate(S t) {
				return asVector(t);
			}
		});
	}

	public static <Q extends Element<Q>, R extends Element<R>> UnivariatePolynomial<Q> minimalPolynomial(R s,
			int degree, Ring<R> extension, Ring<Q> base, MathMap<R, Vector<Q>> asVector) {
		R element = extension.one();
		List<Vector<Q>> generators = new ArrayList<>();
		for (int i = 0; i <= degree; i++) {
			generators.add(asVector.evaluate(element));
			element = extension.multiply(element, s);
			if (generators.size() > 1 && degree % (generators.size() - 1) == 0) {
				MatrixModule<Q> module = new MatrixModule<>(base, degree, generators.size());
				List<Vector<Q>> basis = module.kernelBasis(Matrix.fromColumns(generators));
				if (basis.size() == 1) {
					return base.getUnivariatePolynomialRing()
							.normalize(base.getUnivariatePolynomialRing().getPolynomial(basis.get(0).asList()));
				}
			}
		}
		throw new ArithmeticException("Could not find minimal polynomial for " + s);
	}

	@Override
	public S zero() {
		return zero;
	}

	@Override
	public S one() {
		return one;
	}

	@Override
	public S add(S t1, S t2) {
		return fromSmallDegreePolynomial(polynomials.add(t1.asPolynomial(), t2.asPolynomial()));
	}

	@Override
	public S negative(S t) {
		return fromSmallDegreePolynomial(polynomials.negative(t.asPolynomial()));
	}

	@Override
	public S multiply(S t1, S t2) {
		return fromPolynomial(polynomials.multiply(t1.asPolynomial(), t2.asPolynomial()));
	}

	@Override
	public S inverse(S t) {
		GenericUnivariatePolynomialRing.ExtendedResultantResult<T> eres = polynomials
				.extendedResultant(minimalPolynomial, t.asPolynomial());
		if (eres.getResultant().equals(baseRing.zero())) {
			throw new ArithmeticException("Division by zero: " + t.asPolynomial());
		}
		return fromPolynomial(polynomials.divideScalar(eres.getCoeff2(), eres.getGcd().univariateCoefficient(0)));
	}

	@Override
	public S getRandomElement() {
		return fromPolynomial(polynomials.getRandomElement(degree - 1));
	}

	@Override
	public final Iterator<S> iterator() {
		return new Iterator<>() {
			private Iterator<UnivariatePolynomial<T>> it = polynomials.polynomials(degree() - 1);

			@Override
			public boolean hasNext() {
				return it.hasNext();
			}

			@Override
			public S next() {
				return fromSmallDegreePolynomial(it.next());
			}
		};
	}

	@Override
	public Vector<T> asVector(S s) {
		return polynomials.asVector(s.asPolynomial(), degree - 1);
	}

	@Override
	public S getEmbedding(T t) {
		return fromPolynomial(polynomials.getEmbedding(t));
	}

	@Override
	public S alpha() {
		return alpha;
	}

	public abstract Ext makeExtension(UnivariatePolynomial<T> minimalPolynomial);

	protected abstract Ext asExtensionType();

	private FactorizationResult<Polynomial<T>, T> factorizedMinimalPolynomial() {
		if (factorizedMinimalPolynomial == null) {
			factorizedMinimalPolynomial = baseRing.factorization(minimalPolynomial);
		}
		return factorizedMinimalPolynomial;
	}

	@Override
	public List<Ext> asIrreducibleProduct() {
		if (asIrreducibleProduct == null) {
			FactorizationResult<Polynomial<T>, T> factorizedMipo = factorizedMinimalPolynomial();
			if (factorizedMipo.primeFactors().size() == 1) {
				asIrreducibleProduct = Collections.singletonList(this.asExtensionType());
				return asIrreducibleProduct;
			}
			List<Ext> result = new ArrayList<>();
			for (Polynomial<T> factor : factorizedMipo.primeFactors()) {
				result.add(makeExtension(
						polynomials.toUnivariate(polynomials.power(factor, factorizedMipo.multiplicity(factor)))));
			}
			asIrreducibleProduct = result;
		}
		return asIrreducibleProduct;
	}

	@Override
	public List<S> asIrreducibleProductElement(S s) {
		UnivariatePolynomial<T> asPolynomial = s.asPolynomial();
		List<Ext> asProduct = asIrreducibleProduct();
		List<S> result = new ArrayList<>();
		for (Ext ext : asProduct) {
			result.add(ext.fromPolynomial(asPolynomial));
		}
		return result;
	}

	private List<S> irreducibleProductInverses() {
		if (irreducibleProductInverses == null) {
			List<Ext> exts = asIrreducibleProduct();
			List<S> result = new ArrayList<>();
			for (int i = 0; i < exts.size(); i++) {
				Ext ext = exts.get(i);
				UnivariatePolynomial<T> allButOne = polynomials
						.toUnivariate(polynomials.divideChecked(minimalPolynomial(), ext.minimalPolynomial()));
				result.add(multiply(fromPolynomial(ext.inverse(ext.fromPolynomial(allButOne)).asPolynomial()),
						fromPolynomial(allButOne)));
			}
			irreducibleProductInverses = result;
		}
		return irreducibleProductInverses;
	}

	@Override
	public S fromIrreducibleProductElement(List<S> s) {
		S result = zero();
		List<S> inverses = irreducibleProductInverses();
		for (int i = 0; i < inverses.size(); i++) {
			result = add(result, multiply(fromPolynomial(s.get(i).asPolynomial()), inverses.get(i)));
		}
		return result;
	}

	@Override
	public boolean isUnit(S t) {
		return baseRing.isUnit(polynomials.resultant(t.asPolynomial(), minimalPolynomial));
	}

	@Override
	public boolean isCommutative() {
		return baseRing.isCommutative();
	}

	@Override
	public boolean isIntegral() {
		return baseRing.isIntegral() && polynomials.isIrreducible(minimalPolynomial);
	}

	@Override
	public boolean isReduced() {
		return baseRing.isReduced() && polynomials.squareFreeFactorization(minimalPolynomial).squareFree();
	}

	@Override
	public boolean isIrreducible() {
		return baseRing.isIrreducible() && polynomials.getIdeal(Collections.singletonList(minimalPolynomial))
				.minimalPrimeIdealsOver().size() == 1;
	}

	@Override
	public boolean isZeroDivisor(S t) {
		return baseRing.isZeroDivisor(polynomials.resultant(t.asPolynomial(), minimalPolynomial));
	}

	@Override
	public boolean isEuclidean() {
		return false;
	}

	@Override
	public boolean isUniqueFactorizationDomain() {
		return false;
	}

	@Override
	public FactorizationResult<S, S> uniqueFactorization(S t) {
		throw new ArithmeticException("Not a UFD!");
	}

	@Override
	public boolean isPrincipalIdealDomain() {
		return false;
	}

	@Override
	public boolean isDedekindDomain() {
		return baseRing.isDedekindDomain();
	}

	@Override
	public boolean isDivisible(S dividend, S divisor) {
		return quotientAndRemainder(dividend, divisor).getRemainder().equals(zero());
	}

	@Override
	public QuotientAndRemainderResult<S> quotientAndRemainder(S dividend, S divisor) {
		List<Polynomial<T>> generators = new ArrayList<>();
		generators.add(polynomials.getEmbedding(divisor.asPolynomial()));
		generators.add(polynomials.getEmbedding(minimalPolynomial));
		IdealResult<Polynomial<T>, PolynomialIdeal<T>> divisorIdeal = polynomials.getIdealWithTransforms(generators);
		List<Polynomial<T>> generate = divisorIdeal.getIdeal().generate(dividend.asPolynomial());
		UnivariatePolynomial<T> quotient = polynomials.zero();
		for (int i = 0; i < divisorIdeal.getGeneratorExpressions().size(); i++) {
			quotient = polynomials.add(
					polynomials.multiply(divisorIdeal.getGeneratorExpressions().get(i).get(0), generate.get(i)),
					quotient);
		}
		return new QuotientAndRemainderResult<>(fromPolynomial(quotient),
				fromPolynomial(polynomials.toUnivariate(divisorIdeal.getIdeal().residue(dividend.asPolynomial()))));
	}

	@Override
	public BigInteger euclidMeasure(S t) {
		throw new ArithmeticException("Not Euclidean");
	}

	@Override
	public S projectToUnit(S t) {
		return t;
	}

	@Override
	public Iterable<S> getUnits() {
		throw new UnsupportedOperationException("Not implemented!");
	}

	@Override
	public int krullDimension() {
		return baseRing.krullDimension();
	}

	@Override
	public List<Ideal<S>> maximalPrimeIdealChain(Ideal<S> start) {
		RingExtensionIdeal ideal = (RingExtensionIdeal) start;
		Ideal<T> intersectionIdeal = ideal.polynomialIdeal.intersectToRing();
		List<? extends Ideal<T>> chain = baseRing.maximalPrimeIdealChain(intersectionIdeal);
		List<Ideal<S>> result = new ArrayList<>();
		Ideal<S> prev = start;
		for (Ideal<T> t : chain) {
			Ideal<S> extended = getIdealEmbedding(t, embedding);
			extended = add(extended, prev);
			PrimaryDecompositionResult<S, ? extends Ideal<S>> decomposition = primaryDecomposition(extended);
			prev = decomposition.getRadicals().get(0);
			result.add(prev);
		}
		return result;
	}

	@Override
	public List<Ideal<S>> maximalPrimeIdealChain(Ideal<S> start, Ideal<S> end) {
		RingExtensionIdeal startIdeal = (RingExtensionIdeal) start;
		RingExtensionIdeal endIdeal = (RingExtensionIdeal) end;
		Ideal<T> intersectionIdealStart = startIdeal.polynomialIdeal.intersectToRing();
		Ideal<T> intersectionIdealEnd = endIdeal.polynomialIdeal.intersectToRing();
		List<? extends Ideal<T>> chain = baseRing.maximalPrimeIdealChain(intersectionIdealStart, intersectionIdealEnd);
		List<Ideal<S>> result = new ArrayList<>();
		Ideal<S> prev = start;
		for (Ideal<T> t : chain) {
			Ideal<S> extended = getIdealEmbedding(t, embedding);
			extended = add(extended, prev);
			PrimaryDecompositionResult<S, ? extends Ideal<S>> decomposition = primaryDecomposition(extended);
			for (Ideal<S> radical : decomposition.getRadicals()) {
				if (end.contains(radical)) {
					prev = radical;
					result.add(prev);
					break;
				}
			}
		}
		return result;
	}

	@Override
	public CoordinateRing<T> asCoordinateRing() {
		return polynomials.getIdeal(Collections.singletonList(minimalPolynomial)).divideOut();
	}

	@Override
	public PolynomialRingAsCoordinateRing<T, S> asCoordinateRing(PolynomialRing<S> polynomialRing) {
		PolynomialRing<T> polynomialsForCoordinateRing = AbstractPolynomialRing.getPolynomialRing(baseRing,
				polynomialRing.numberOfVariables() + 1, polynomialRing.getComparator());
		Polynomial<T> embeddedMinimalPolynomial = polynomialsForCoordinateRing.getEmbedding(minimalPolynomial,
				new int[] { polynomialRing.numberOfVariables() });
		CoordinateRing<T> coordinateRing = polynomialsForCoordinateRing
				.getIdeal(Collections.singletonList(embeddedMinimalPolynomial)).divideOut();
		return new PolynomialRingAsCoordinateRing<>(polynomialRing, coordinateRing, new MathMap<>() {

			@Override
			public CoordinateRingElement<T> evaluate(Polynomial<S> t) {
				Map<Monomial, T> coefficients = new TreeMap<>();
				for (Monomial m : t.monomials()) {
					S coefficient = t.coefficient(m);
					UnivariatePolynomial<T> asPolynomial = coefficient.asPolynomial();
					int[] exponents = Arrays.copyOf(m.exponents(), polynomialRing.numberOfVariables() + 1);
					for (int i = 0; i <= asPolynomial.degree(); i++) {
						exponents[polynomialRing.numberOfVariables()] = i;
						Monomial m2 = polynomialsForCoordinateRing.getMonomial(exponents);
						coefficients.put(m2, asPolynomial.univariateCoefficient(i));
					}
				}
				return coordinateRing.getEmbedding(polynomialsForCoordinateRing.getPolynomial(coefficients));
			}
		}, new MathMap<>() {

			@Override
			public Polynomial<S> evaluate(CoordinateRingElement<T> t) {
				Polynomial<T> asPolynomial = t.getElement();
				Map<Monomial, S> coefficients = new TreeMap<>();
				for (Monomial m : asPolynomial.monomials()) {
					int[] exponents = Arrays.copyOf(m.exponents(), polynomialRing.numberOfVariables());
					Monomial m2 = polynomialRing.getMonomial(exponents);
					S coefficient = coefficients.getOrDefault(m2, zero());
					coefficients.put(m2, add(coefficient, scalarMultiply(asPolynomial.coefficient(m),
							power(alpha(), m.exponents()[polynomialRing.numberOfVariables()]))));
				}
				return polynomialRing.getPolynomial(coefficients);
			}
		});
	}

	@Override
	public PolynomialRingAsCoordinateRing<T, S> asCoordinateRing(int numberOfVariables) {
		return asCoordinateRing(numberOfVariables, Monomial.GREVLEX);
	}

	@Override
	public PolynomialRingAsCoordinateRing<T, S> asCoordinateRing(int numberOfVariables,
			Comparator<Monomial> comparator) {
		return asCoordinateRing(AbstractPolynomialRing.getPolynomialRing(this, numberOfVariables, comparator));
	}

	@Override
	public ExtensionCoordinateRing<T, S> asCoordinateRing(CoordinateRing<S> coordinateRing) {
		PolynomialRingAsCoordinateRing<T, S> polynomialRing = asCoordinateRing(coordinateRing.getPolynomialRing());
		List<CoordinateRingElement<T>> generators = new ArrayList<>();
		for (Polynomial<S> generator : coordinateRing.getIdeal().generators()) {
			generators.add(polynomialRing.getIsomorphism().evaluate(generator));
		}
		CoordinateIdeal<T> ideal = polynomialRing.getCoordinateRing().getIdeal(generators);
		CoordinateRing<T> baseCoordinateRing = ideal.divideOut();
		return new ExtensionCoordinateRing<>(coordinateRing, baseCoordinateRing, new MathMap<>() {

			@Override
			public CoordinateRingElement<T> evaluate(CoordinateRingElement<S> t) {
				return baseCoordinateRing
						.getEmbedding(polynomialRing.getIsomorphism().evaluate(t.getElement()).getElement());
			}
		}, new MathMap<>() {

			@Override
			public CoordinateRingElement<S> evaluate(CoordinateRingElement<T> t) {
				return coordinateRing.getEmbedding(polynomialRing.getInverseIsomorphism()
						.evaluate(polynomialRing.getCoordinateRing().getEmbedding(t.getElement())));
			}
		});
	}

	@Override
	public IdealResult<S, RingExtensionIdeal> getIdealWithTransforms(List<S> generators) {
		List<Polynomial<T>> asPolynomials = new ArrayList<>();
		asPolynomials.add(minimalPolynomial());
		for (S s : generators) {
			asPolynomials.add(s.asPolynomial());
		}
		IdealResult<Polynomial<T>, PolynomialIdeal<T>> polynomialIdealResult = polynomials
				.getIdealWithTransforms(asPolynomials);
		RingExtensionIdeal ideal = new RingExtensionIdeal(polynomialIdealResult.getIdeal());
		List<List<S>> expressions = new ArrayList<>();
		for (List<Polynomial<T>> expression : polynomialIdealResult.getGeneratorExpressions()) {
			List<S> fromPolynomial = new ArrayList<>();
			for (Polynomial<T> p : expression) {
				fromPolynomial.add(fromPolynomial(polynomials.toUnivariate(p)));
			}
			expressions.add(fromPolynomial);
		}
		List<Vector<S>> syzygies = new ArrayList<>();
		for (Vector<Polynomial<T>> syzygy : polynomialIdealResult.getSyzygies()) {
			List<S> fromPolynomial = new ArrayList<>();
			for (Polynomial<T> p : syzygy.asList()) {
				fromPolynomial.add(fromPolynomial(polynomials.toUnivariate(p)));
			}
			syzygies.add(new Vector<>(fromPolynomial));
		}
		return new IdealResult<>(expressions, generators, ideal, syzygies);
	}

	@Override
	public RingExtensionIdeal intersect(Ideal<S> t1, Ideal<S> t2) {
		RingExtensionIdeal s1 = (RingExtensionIdeal) t1;
		RingExtensionIdeal s2 = (RingExtensionIdeal) t2;
		return new RingExtensionIdeal(polynomials.intersect(s1.polynomialIdeal, s2.polynomialIdeal));
	}

	@Override
	public RingExtensionIdeal radical(Ideal<S> t) {
		RingExtensionIdeal s = (RingExtensionIdeal) t;
		return new RingExtensionIdeal(polynomials.radical(s.polynomialIdeal));
	}

	@Override
	public Ideal<S> getNilRadical() {
		FactorizationResult<Polynomial<T>, T> factorizedMipo = factorizedMinimalPolynomial();
		UnivariatePolynomial<T> nilRadical = polynomials.one();
		for (Polynomial<T> factor : factorizedMipo.primeFactors()) {
			nilRadical = polynomials.multiply(factor, nilRadical);
		}
		List<S> generators = new ArrayList<>();
		generators.add(fromPolynomial(nilRadical));
		for (T generator : baseRing.getNilRadical().generators()) {
			generators.add(getEmbedding(generator));
		}
		return getIdeal(generators);
	}

	private class RingExtensionIdeal extends AbstractIdeal<S> implements Ideal<S> {
		private PolynomialIdeal<T> polynomialIdeal;
		private List<S> generators;

		private RingExtensionIdeal(PolynomialIdeal<T> polynomialIdeal) {
			super(AbstractAlgebraicRingExtension.this);
			this.polynomialIdeal = polynomialIdeal;
			this.generators = new ArrayList<>();
			for (Polynomial<T> g : polynomialIdeal.generators()) {
				this.generators.add(fromPolynomial(polynomials.toUnivariate(g)));
			}
		}

		@Override
		public boolean isFinite() {
			return baseRing.isFinite()
					|| polynomialIdeal.equals(polynomials.getIdeal(Collections.singletonList(minimalPolynomial())));
		}

		@Override
		public BigInteger getNumberOfElements() throws InfinityException {
			throw new UnsupportedOperationException("Not implemented!");
		}

		@Override
		public boolean isPrimary() {
			return polynomialIdeal.isPrimary();
		}

		@Override
		public boolean isPrime() {
			return polynomialIdeal.isPrime();
		}

		@Override
		public boolean isMaximal() {
			return polynomialIdeal.isMaximal();
		}

//		@Override
//		public List<Vector<S>> nonTrivialCombinations(List<S> s) {
//			List<Polynomial<T>> asPolynomials = new ArrayList<>();
//			for (S e : s) {
//				asPolynomials.add(asPolynomial(e));
//			}
//			List<Vector<S>> result = new ArrayList<>();
//			List<Vector<Polynomial<T>>> combinations = polynomials.syzygyProblem(asPolynomials);
//			for (Vector<Polynomial<T>> combination : combinations) {
//				List<S> row = new ArrayList<>();
//				for (Polynomial<T> e : combination.asList()) {
//					row.add(fromPolynomial(polynomials.toUnivariate(e)));
//				}
//				result.add(new Vector<>(row));
//			}
//			return result;
//		}

		@Override
		public List<S> generators() {
			return generators;
		}

		@Override
		public List<S> generate(S t) {
			List<Polynomial<T>> generated = polynomialIdeal.generate(t.asPolynomial());
			List<S> result = new ArrayList<>();
			for (Polynomial<T> g : generated) {
				result.add(fromPolynomial(polynomials.toUnivariate(g)));
			}
			return result;
		}

		@Override
		public S residue(S t) {
			return fromPolynomial(polynomials.toUnivariate(polynomialIdeal.residue(t.asPolynomial())));
		}

		@Override
		public boolean contains(S t) {
			return polynomialIdeal.contains(t.asPolynomial());
		}

		@Override
		public Value maximumPowerContains(S t) {
			return polynomialIdeal.maximumPowerContains(t.asPolynomial());
		}

	}
}
