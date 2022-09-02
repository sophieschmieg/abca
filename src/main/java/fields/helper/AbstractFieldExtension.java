package fields.helper;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;

import fields.exceptions.InfinityException;
import fields.helper.GenericAlgebraicRingExtension.GenericAlgebraicExtensionElement;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals.Fraction;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.BilinearMap;
import fields.interfaces.DedekindRing;
import fields.interfaces.DiscreteValuationField;
import fields.interfaces.DiscreteValuationRing;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.FieldExtension;
import fields.interfaces.GlobalField;
import fields.interfaces.Ideal;
import fields.interfaces.LocalRing;
import fields.interfaces.MathMap;
import fields.interfaces.Module;
import fields.interfaces.Polynomial;
import fields.interfaces.Ring;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.LocalRingImplementation;
import fields.local.TrivialDiscreteValuationField;
import fields.local.Value;
import fields.polynomials.FieldExtensionUnivariatePolynomialRing;
import fields.vectors.FiniteVectorSpace;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;
import fields.vectors.pivot.PivotStrategy;
import fields.vectors.pivot.TrivialPivotStrategy;
import util.ConstantMap;
import util.Identity;
import util.Pair;
import util.SingletonSortedMap;

public abstract class AbstractFieldExtension<T extends Element<T>, S extends AlgebraicExtensionElement<T, S>, Ext extends FieldExtension<T, S, Ext>>
		extends AbstractAlgebraicRingExtension<T, S, Ext> implements FieldExtension<T, S, Ext>, LocalRing<S, S, S> {
	private Field<T> baseField;
	private UnivariatePolynomialRing<T> polynomials;
	private UnivariatePolynomial<T> minimalPolynomial;
	private int degree;
	private FiniteVectorSpace<T> asVectorSpace;
	private Polynomial<T> genericNorm;
	private FieldExtensionUnivariatePolynomialRing<T, S, Ext> univariatePolynomialRing;

	public AbstractFieldExtension(UnivariatePolynomial<T> minimalPolynomial, Field<T> baseField, String variableName) {
		super(minimalPolynomial, baseField, variableName);
		this.baseField = baseField;
		this.polynomials = baseField.getUnivariatePolynomialRing();
		this.minimalPolynomial = minimalPolynomial;
		this.degree = this.minimalPolynomial.degree();
		this.asVectorSpace = new FiniteVectorSpace<>(this.baseField, degree);
	}

	public AbstractFieldExtension(Field<T> baseField) {
		this(baseField.getUnivariatePolynomialRing().getVar(), baseField, "x");
	}

	public Polynomial<T> genericNorm() {
		if (genericNorm == null) {
			genericNorm = genericFreeModule().matrixAlgebra().determinant(genericMatrix());
		}
		return genericNorm;
	}

	@Override
	public String toString() {
		if (degree() == 1) {
			return baseField.toString();
		}
		return polynomials.toString() + "/(" + minimalPolynomial.toString() + ")";
	}

	@Override
	public final FiniteVectorSpace<T> asVectorSpace() {
		return asVectorSpace;
	}

	public final S getUnitVector(int index) {
		return power(alpha(), index - 1);
	}

	@Override
	public final List<S> getBasis() {
		return getModuleGenerators();
	}

	protected abstract S fromSmallDegreePolynomial(UnivariatePolynomial<T> polynomial);

	@Override
	public final Field<T> getField() {
		return baseField;
	}

	@Override
	public Field<T> getBaseField() {
		return baseField;
	}

	@Override
	public final int degreeOver(FieldEmbedding<T, S, Ext> base) {
		if (!base.getField().equals(this)) {
			throw new ArithmeticException("Not embedded here!");
		}
		return degree() / base.getEmbeddedField().degree();
	}

	@Override
	public final boolean isFiniteExtension() {
		return true;
	}

	@Override
	public final boolean isSeparable() {
		return !polynomials.derivative(minimalPolynomial()).equals(polynomials.zero());
	}

	@Override
	public final boolean isNormal() {
		return roots(getUnivariatePolynomialRing().getEmbedding(minimalPolynomial(), getEmbeddingMap()))
				.size() == degree();
	}

	@Override
	public final boolean isGalois() {
		return isFiniteExtension() && isSeparable() && isNormal();
	}

	@Override
	public boolean isCyclic() {
		return false;
	}

	@Override
	public FieldAutomorphism<T, S, Ext> cyclicGenerator() {
		throw new UnsupportedOperationException("Not implemented");
	}

	@Override
	public int dimension() {
		return degree();
	}

	@Override
	public DedekindRing<S, S, S> asDedekindRing() {
		return this;
	}

	@Override
	public S divide(S dividend, S divisor) {
		return multiply(dividend, inverse(divisor));
	}

	@Override
	public final Matrix<S> asMatrixOver(S t, FieldEmbedding<T, S, Ext> base) {
		return asMatrix(asVectorOver(t, base), minimalPolynomialOver(base), base.getEmbeddedField().zero(),
				base.asVectorSpace());
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
	public T norm(S t) {
		return matrixAlgebra().determinant(asMatrix(t));
	}

	@Override
	public final T trace(S t) {
		return matrixAlgebra().trace(asMatrix(t));
	}

	@Override
	public S normOver(S t, FieldEmbedding<T, S, Ext> base) {
		if (!base.getField().equals(this)) {
			throw new ArithmeticException("Not embedded here!");
		}
		return base.norm(t);
	}

	@Override
	public S traceOver(S t, FieldEmbedding<T, S, Ext> base) {
		if (!base.getField().equals(this)) {
			throw new ArithmeticException("Not embedded here!");
		}
		return base.trace(t);
	}

	@Override
	public T traceForm(S s1, S s2) {
		return trace(multiply(s1, s2));
	}

	@Override
	public BilinearMap<S, T> traceForm() {
		return new BilinearMap<>() {
			@Override
			public T evaluate(S t1, S t2) {
				return traceForm(t1, t2);
			}
		};
	}

	@Override
	public S traceFormOver(S s1, S s2, FieldEmbedding<T, S, Ext> base) {
		return traceOver(multiply(s1, s2), base);
	}

	@Override
	public BilinearMap<S, S> traceFormOver(FieldEmbedding<T, S, Ext> base) {
		return new BilinearMap<>() {
			@Override
			public S evaluate(S t1, S t2) {
				return traceFormOver(t1, t2, base);
			}
		};
	}

	@Override
	public Matrix<T> traceFormMatrix() {
		return Matrix.fromBilinearMap(this, traceForm());
	}

	@Override
	public Matrix<S> traceFormMatrixOver(FieldEmbedding<T, S, Ext> base) {
		return Matrix.fromBilinearMap(base.asVectorSpace(), new BilinearMap<>() {
			@Override
			public S evaluate(Vector<S> t1, Vector<S> t2) {
				return traceFormOver(base).evaluate(base.fromVector(t1), base.fromVector(t2));
			}
		});
	}

	@Override
	public final UnivariatePolynomial<S> minimalPolynomialOver(FieldEmbedding<T, S, Ext> base) {
		if (!base.getField().equals(this)) {
			throw new ArithmeticException("Not embedded here");
		}
		return minimalPolynomial(alpha(), degreeOver(base), this, base.getEmbeddedField(), new MathMap<>() {
			@Override
			public Vector<S> evaluate(S t) {
				return asVectorOver(t, base);
			}
		});
	}

	@Override
	public final UnivariatePolynomial<S> minimalPolynomialOver(S s, FieldEmbedding<T, S, Ext> base) {
		if (!base.getField().equals(this)) {
			throw new ArithmeticException("Not embedded here");
		}
		return minimalPolynomial(s, degreeOver(base), this, base.getEmbeddedField(), new MathMap<>() {
			@Override
			public Vector<S> evaluate(S t) {
				return asVectorOver(t, base);
			}
		});
	}

	@Override
	public Vector<S> asVectorOver(S s, FieldEmbedding<T, S, Ext> base) {
		if (!base.getField().equals(this)) {
			throw new ArithmeticException("Not embedded here");
		}
		return base.asVector(s);
	}

	public static <Q extends Element<Q>, R extends Element<R>> UnivariatePolynomial<Q> minimalPolynomial(R s,
			int degree, Ring<R> extension, Field<Q> base, MathMap<R, Vector<Q>> asVector) {
		R element = extension.one();
		List<Vector<Q>> generators = new ArrayList<>();
		for (int i = 0; i <= degree; i++) {
			generators.add(asVector.evaluate(element));
			if (generators.size() > 1 && degree % (generators.size() - 1) == 0) {
				MatrixModule<Q> module = new MatrixModule<>(base, generators.get(0).dimension(), generators.size());
				List<Vector<Q>> basis = module.kernelBasis(Matrix.fromColumns(generators));
				if (basis.size() == 1) {
					return base.getUnivariatePolynomialRing()
							.normalize(base.getUnivariatePolynomialRing().getPolynomial(basis.get(0).asList()));
				}
			}
			element = extension.multiply(element, s);
		}
		throw new ArithmeticException("Could not find minimal polynomial for " + s);
	}

	@Override
	public List<S> conjugates(S s) {
		List<S> conjugates = new ArrayList<>();
		conjugates.addAll(
				roots(getUnivariatePolynomialRing().getEmbedding(minimalPolynomial(s), getEmbeddingMap())).keySet());
		return conjugates;
	}

	@Override
	public GaloisGroup<T, S, Ext> galoisGroup() {
		throw new UnsupportedOperationException("Not implemented");
	}

	@Override
	public S hilbert90(S s) {
		if (!isCyclic() || !norm(s).equals(getBaseField().one())) {
			throw new ArithmeticException("Hilbert 90 preconditions not met");
		}
		Matrix<T> sigma = Matrix.fromEndomorphism(this, cyclicGenerator());
		Matrix<T> asMatrix = asMatrix(s);
		MatrixAlgebra<T> algebra = matrixAlgebra();
		return fromVector(algebra.kernelBasis(algebra.subtract(sigma, asMatrix)).get(0));
	}

	private MathMap<GenericAlgebraicExtensionElement<S>, Vector<T>> asVectorAlphaBetaGeneric(int degree) {
		return new MathMap<>() {
			@Override
			public Vector<T> evaluate(GenericAlgebraicExtensionElement<S> t) {
				List<T> asVector = new ArrayList<>();
				for (int i = 0; i < degree; i++) {
					asVector.addAll(
							AbstractFieldExtension.this.asVector(t.asPolynomial().univariateCoefficient(i)).asList());
				}
				return new Vector<>(asVector);
			}
		};
	}

	protected UnivariatePolynomial<T> minimalPolynomialOverPrimeGeneric(int degree, GenericAlgebraicRingExtension<S> cr,
			GenericAlgebraicExtensionElement<S> s) {
		return minimalPolynomial(s, degree * degree(), cr, baseField, asVectorAlphaBetaGeneric(degree));
	}

	@Override
	public abstract Ext makeExtension(UnivariatePolynomial<T> minimalPolynomial);

	protected abstract Ext asExtensionType();

	@Override
	public FieldEmbedding<T, S, Ext> getEmbeddedExtension(UnivariatePolynomial<S> minimalPolynomial) {
		if (getUnivariatePolynomialRing().derivative(minimalPolynomial).equals(getUnivariatePolynomialRing().zero())) {
			throw new ArithmeticException("Inseparable field extensions are not supported at the moment");
		}
		if (minimalPolynomial.degree() == 1) {
			minimalPolynomial = getUnivariatePolynomialRing().normalize(minimalPolynomial);
			return new FieldEmbedding<>(asExtensionType(), asExtensionType(), alpha(),
					negative(minimalPolynomial.univariateCoefficient(0)));
		}
		GenericAlgebraicRingExtension<S> cr = new GenericAlgebraicRingExtension<>(minimalPolynomial, this);
		int degree = minimalPolynomial.degree();
		GenericAlgebraicExtensionElement<S> alpha = cr.alpha();
		GenericAlgebraicExtensionElement<S> beta = cr.getEmbedding(alpha());
		GenericAlgebraicExtensionElement<S> gamma = null;
		Ext extension = null;
		if (isFinite()) {
			for (GenericAlgebraicExtensionElement<S> element : cr) {
				gamma = element;
				UnivariatePolynomial<T> minimalPolynomialOverPrime = minimalPolynomialOverPrimeGeneric(degree, cr,
						gamma);
				if (minimalPolynomialOverPrime.degree() == degree * degree()) {
					extension = makeExtension(minimalPolynomialOverPrime);
					break;
				}
			}
		} else {
			BigInteger i = BigInteger.ZERO;
			boolean randomize = false;
			while (true) {
				GenericAlgebraicExtensionElement<S> lambda;
				if (randomize) {
					lambda = cr.getEmbedding(getRandomElement());
				} else {
					lambda = cr.getInteger(i);
				}
				gamma = cr.add(alpha, cr.multiply(lambda, beta));
				UnivariatePolynomial<T> minimalPolynomialOverPrime = minimalPolynomialOverPrimeGeneric(degree, cr,
						gamma);
				if (minimalPolynomialOverPrime.degree() == degree * degree()) {
					extension = makeExtension(minimalPolynomialOverPrime);
					break;
				}
				if (characteristic().equals(BigInteger.ZERO) || characteristic().compareTo(i.add(BigInteger.ONE)) > 0) {
					i = i.add(BigInteger.ONE);
				} else {
					randomize = true;
				}
			}
		}
		MathMap<GenericAlgebraicExtensionElement<S>, Vector<T>> asVector = asVectorAlphaBetaGeneric(degree);
		List<Vector<T>> gammaBasis = new ArrayList<>();
		GenericAlgebraicExtensionElement<S> element = cr.getEmbedding(one());
		for (int i = 0; i < degree * degree(); i++) {
			gammaBasis.add(asVector.evaluate(element));
			element = cr.multiply(element, gamma);
		}
		Matrix<T> gammaBaseToAlphaBetaBase = Matrix.fromColumns(gammaBasis);
		Matrix<T> alphaBetaBaseToGammaBase = extension.matrixAlgebra().inverse(gammaBaseToAlphaBetaBase);
		return new FieldEmbedding<T, S, Ext>(extension, this.asExtensionType(),
				extension.fromVector(
						extension.matrixAlgebra().multiply(alphaBetaBaseToGammaBase, asVector.evaluate(beta))),
				extension.fromVector(
						extension.matrixAlgebra().multiply(alphaBetaBaseToGammaBase, asVector.evaluate(alpha))));
	}

	@Override
	public Extension<S, T, S, Ext> getExtension(UnivariatePolynomial<S> minimalPolynomial) {
		return getEmbeddedExtension(minimalPolynomial).asExtension();
	}

	@Override
	public SplittingFieldResult<T, S, Ext> getSplittingField(UnivariatePolynomial<S> minimalPolynomial) {
		if (minimalPolynomial.degree() == 0) {
			return new SplittingFieldResult<T, S, Ext>(new FieldEmbedding<>(asExtensionType()), Collections.emptyList(),
					minimalPolynomial);
		}
		UnivariatePolynomialRing<S> polynomialRing = getUnivariatePolynomialRing();
		FactorizationResult<Polynomial<S>, S> factorization = factorization(minimalPolynomial);
		UnivariatePolynomial<S> firstFactor = polynomialRing.toUnivariate(factorization.firstPrimeFactor());
		FieldEmbedding<T, S, Ext> firstEmbedding;
		S firstConjugate;
		if (firstFactor.degree() == 1) {
			firstEmbedding = new FieldEmbedding<>(asExtensionType());
			firstConjugate = negative(firstFactor.univariateCoefficient(0));
		} else {
			firstEmbedding = getEmbeddedExtension(firstFactor);
			firstConjugate = firstEmbedding.getGenerator();
		}
		polynomialRing = firstEmbedding.getField().getUnivariatePolynomialRing();
		UnivariatePolynomial<S> firstLinear = polynomialRing.toUnivariate(
				polynomialRing.subtract(polynomialRing.getVar(), polynomialRing.getEmbedding(firstConjugate)));
		UnivariatePolynomial<S> remaining = polynomialRing
				.toUnivariate(polynomialRing.divideChecked(polynomialRing.getEmbedding( minimalPolynomial, firstEmbedding.getEmbeddingMap()), firstLinear));
		SplittingFieldResult<T, S, Ext> remainingSplittingField = firstEmbedding.getField().getSplittingField(remaining);
		List<S> conjugates = new ArrayList<>();
		conjugates.add(remainingSplittingField.getExtension().getEmbedding(firstConjugate));
		conjugates.addAll(remainingSplittingField.getConjugates());
		return new SplittingFieldResult<>(new FieldEmbedding<>(firstEmbedding, remainingSplittingField.getExtension()), conjugates, minimalPolynomial);
	}

	@Override
	public boolean isIrreducible(UnivariatePolynomial<S> t) {
		return factorization(t).isIrreducible();
	}

	@Override
	public boolean isSubModuleMember(MatrixModule<S> module, Matrix<S> m, Vector<S> b) {
		MatrixModule<S>.LDUPResult ldup = module.ldup(m);
		// LD^-1Ux=P^-1b
		Vector<S> rhs = module.permuteVector(ldup.getInversePermutation(), b);
		// D^-1Ux=L^-1P^-1b
		List<S> solution = new ArrayList<>();
		for (int i = 0; i < m.rows(); i++) {
			S adjustedRhs = rhs.get(i + 1);
			for (int k = 0; k < i; k++) {
				adjustedRhs = subtract(adjustedRhs,
						multiply(ldup.getLowerTriangle().entry(i + 1, k + 1), solution.get(k)));
			}
			solution.add(divide(adjustedRhs, ldup.getLowerTriangle().entry(i + 1, i + 1)));
		}
		for (int i = ldup.getRank(); i < b.dimension(); i++) {
			if (!solution.get(i).equals(zero())) {
				return false;
			}
		}
		return true;
	}

	@Override
	public Vector<S> asSubModuleMember(MatrixModule<S> module, Matrix<S> m, Vector<S> b) {
		MatrixModule<S>.LDUPResult ldup = module.ldup(m);
		// PLD^-1Ux=Ax=b
		// LD^-1Ux=P^-1b
		Vector<S> rhs = module.permuteVector(ldup.getInversePermutation(), b);
		// D^-1Ux=L^-1P^-1b
		List<S> solution = new ArrayList<>();
		for (int i = 0; i < m.rows(); i++) {
			S adjustedRhs = rhs.get(i + 1);
			for (int k = 0; k < i; k++) {
				adjustedRhs = subtract(adjustedRhs,
						multiply(ldup.getLowerTriangle().entry(i + 1, k + 1), solution.get(k)));
			}
			solution.add(divide(adjustedRhs, ldup.getLowerTriangle().entry(i + 1, i + 1)));
		}
		// Ux = DL^-1P^-1b
		List<S> adjustedSolution = new ArrayList<>();
		for (int i = 0; i < m.rows(); i++) {
			adjustedSolution.add(multiply(ldup.getInverseDiagonal().entry(i + 1, i + 1), solution.get(i)));
		}
		// x = U^-1DL^-1P^-1b
		List<S> reverseSolution = new ArrayList<>();
		int j = m.columns();
		for (int i = m.rows(); i > 0; i--) {
			if (!ldup.getSteps().containsKey(i)) {
				continue;
			}
			while (ldup.getSlips().contains(j)) {
				reverseSolution.add(zero());
				j--;
			}
			S adjustedRhs = adjustedSolution.get(i - 1);
			for (int k = 0; k < reverseSolution.size(); k++) {
				adjustedRhs = subtract(adjustedRhs,
						multiply(reverseSolution.get(k), ldup.getUpperTriangle().entry(i, m.columns() - k)));
			}
			reverseSolution.add(divide(adjustedRhs, ldup.getUpperTriangle().entry(i, j)));
			j--;
		}
		while (ldup.getSlips().contains(j)) {
			reverseSolution.add(zero());
			j--;
		}
		Collections.reverse(reverseSolution);
		return new Vector<>(reverseSolution);
	}

	@Override
	public List<Vector<S>> syzygyProblem(MatrixModule<S> module, Matrix<S> m) {
		MatrixModule<S>.LDUPResult ldup = module.ldup(m);
		List<Vector<S>> result = new ArrayList<>();
		for (int slip : ldup.getSlips()) {
			List<S> reverseSolution = new ArrayList<>();
			for (int i = 0; i < m.columns() - slip; i++) {
				reverseSolution.add(zero());
			}
			reverseSolution.add(one());
			int index = slip - 1;
			int firstStepRow = -1;
			for (int i = ldup.getRank(); i > 0; i--) {
				int step = ldup.getSteps().get(i);
				if (step > slip) {
					continue;
				}
				if (firstStepRow < 0) {
					firstStepRow = i;
				}
				S slipValue = ldup.getUpperTriangle().entry(i, slip);
				for (; index > step; index--) {
					reverseSolution.add(zero());
				}
				for (int j = firstStepRow; j > i; j--) {
					int jStep = ldup.getSteps().get(j);
					slipValue = add(
							multiply(ldup.getUpperTriangle().entry(i, jStep), reverseSolution.get(m.columns() - jStep)),
							slipValue);
				}
				index--;
				reverseSolution.add(divide(negative(slipValue), ldup.getUpperTriangle().entry(i, step)));
			}
			for (; index > 0; index--) {
				reverseSolution.add(zero());
			}
			Collections.reverse(reverseSolution);
			result.add(new Vector<>(reverseSolution));
		}
		return result;
	}

	@Override
	public List<Vector<S>> simplifySubModuleGenerators(MatrixModule<S> module, Matrix<S> m) {
		List<Vector<S>> asVectors = new ArrayList<>();
		MatrixModule<S>.LDUPResult gauss = module.ldup(m);
		for (int j = 0; j < m.columns(); j++) {
			if (gauss.getSlips().contains(j + 1)) {
				continue;
			}
			asVectors.add(m.column(j + 1));
		}
		return asVectors;
	}

	@Override
	public int krullDimension() {
		return 0;
	}

	@Override
	public List<Ideal<S>> maximalPrimeIdealChain() {
		return Collections.singletonList(getZeroIdeal());
	}

	@Override
	public List<Ideal<S>> maximalPrimeIdealChain(Ideal<S> start) {
		return Collections.singletonList(getZeroIdeal());
	}

	@Override
	public List<Ideal<S>> maximalPrimeIdealChain(Ideal<S> start, Ideal<S> end) {
		return Collections.singletonList(getZeroIdeal());
	}

	@Override
	public boolean isUniqueFactorizationDomain() {
		return true;
	}

	@Override
	public FactorizationResult<S, S> uniqueFactorization(S t) {
		return new FactorizationResult<>(t, Collections.emptySortedMap());
	}

	@Override
	public boolean isIrreducible(S t) {
		throw new ArithmeticException("not a non zero non unit! (it's a field!)");
	}

	@Override
	public boolean isPrime(S t) {
		throw new ArithmeticException("not a non zero non unit! (it's a field!)");
	}

	@Override
	public List<S> factors(S t) {
		return Collections.singletonList(one());
	}

	@Override
	public List<S> adicDevelopment(S t, S base) {
		throw new ArithmeticException("No adic development over fields!");
	}

	@Override
	public boolean isPrincipalIdealDomain() {
		return true;
	}

	@Override
	public boolean isDedekindDomain() {
		return true;
	}

	@Override
	public FactorizationResult<Ideal<S>, Ideal<S>> idealFactorization(Ideal<S> t) {
		if (t.equals(getZeroIdeal())) {
			return new FactorizationResult<>(getUnitIdeal(), SingletonSortedMap.map(getZeroIdeal(), 1));
		}
		return new FactorizationResult<>(getUnitIdeal(), new TreeMap<>());
	}

	@Override
	public PrimaryDecompositionResult<S, Ideal<S>> primaryDecomposition(Ideal<S> t) {
		if (t.equals(getZeroIdeal())) {
			return new PrimaryDecompositionResult<>(Collections.singletonList(getZeroIdeal()),
					Collections.singletonList(getZeroIdeal()));
		}
		return new PrimaryDecompositionResult<>(Collections.emptyList(), Collections.emptyList());
	}

	@Override
	public ModuloMaximalIdealResult<S, S, Field<S>, Ideal<S>, Field<S>> moduloMaximalIdeal(Ideal<S> ideal) {
		if (!ideal.equals(getZeroIdeal())) {
			throw new ArithmeticException("Not a maximal ideal!");
		}
		return new ModuloMaximalIdealResult<>(this, ideal, this, new Identity<>(), new Identity<>());
	}

	@Override
	public ModuloIdealResult<S, S> moduloIdeal(Ideal<S> ideal) {
		if (ideal.contains(one())) {
			throw new ArithmeticException("Not a proper ideal!");
		}
		return new ModuloIdealResult<>(this, ideal, this, new Identity<>(), new Identity<>());
	}

	@Override
	public <U extends Element<U>> Ideal<S> getIdealEmbedding(Ideal<U> t, MathMap<U, S> map) {
		for (U generator : t.generators()) {
			if (!map.evaluate(generator).equals(zero())) {
				return getUnitIdeal();
			}
		}
		return getZeroIdeal();
	}

	@Override
	public Ideal<S> getNilRadical() {
		return getZeroIdeal();
	}

//	@Override
//	public IdealResult<S, FieldIdeal<S>> getIdealWithTransforms(List<S> generators) {
//		FieldIdeal<S> ideal = new FieldIdeal<>(generators, this);
//		List<S> expression = new ArrayList<>();
//		boolean found = false;
//		for (S g : generators) {
//			if (!found && !g.equals(zero())) {
//				found = true;
//				expression.add(inverse(g));
//			} else {
//				expression.add(zero());
//			}
//		}
//		return new IdealResult<>(Collections.singletonList(expression), generators, ideal);
//	}
//
//	public Ideal<S> add(Ideal<S> t1, Ideal<S> t2) {
//		return getIdeal(Collections.singletonList(
//				t1.generators().get(0).equals(zero()) && t2.generators().get(0).equals(zero()) ? zero() : one()));
//	}
//
//	public Ideal<S> multiply(Ideal<S> t1, Ideal<S> t2) {
//		return getIdeal(Collections.singletonList(
//				t1.generators().get(0).equals(zero()) || t2.generators().get(0).equals(zero()) ? zero() : one()));
//	}
//
//	public Ideal<S> intersect(Ideal<S> t1, Ideal<S> t2) {
//		return getIdeal(Collections.singletonList(
//				t1.generators().get(0).equals(zero()) || t2.generators().get(0).equals(zero()) ? zero() : one()));
//	}
//
//	@Override
//	public Ideal<S> radical(Ideal<S> t) {
//		return t;
//	}

	@Override
	public S getFraction(Fraction t) {
		return divide(getInteger(t.getNumerator()), getInteger(t.getDenominator()));
	}

	public Iterable<S> getNonZeroElements() throws InfinityException {
		return new Iterable<S>() {
			private Iterable<S> elements = AbstractFieldExtension.this;

			@Override
			public Iterator<S> iterator() {
				return new Iterator<S>() {
					private Iterator<S> it = elements.iterator();

					@Override
					public boolean hasNext() {
						return it.hasNext();
					}

					@Override
					public S next() {
						S result = null;
						do {
							if (it.hasNext())
								result = it.next();
						} while (result.equals(AbstractFieldExtension.this.zero()));
						return result;
					}
				};
			}

		};
	}

	@Override
	public S primitiveRoot() {
		throw new UnsupportedOperationException();
	}

	@Override
	public BigInteger getNumberOfUnits() {
		return this.getNumberOfElements().subtract(BigInteger.ONE);
	}

	@Override
	public boolean isUnit(S t) {
		return !t.equals(this.zero());
	}

	@Override
	public S projectToUnit(S t) {
		if (t.equals(zero())) {
			return one();
		}
		return t;
	}

	@Override
	public S upToUnit(S t) {
		if (t.equals(zero())) {
			return zero();
		}
		return one();
	}

	@Override
	public boolean isIntegral() {
		return true;
	}

	@Override
	public boolean isReduced() {
		return true;
	}

	@Override
	public boolean isIrreducible() {
		return true;
	}

	@Override
	public Value valuation(S t, Ideal<S> maximalIdeal) {
		return t.equals(zero()) ? Value.INFINITY : Value.ZERO;
	}

	@Override
	public FieldOfFractionsResult<S, S> fieldOfFractions() {
		return new FieldOfFractionsResult<>(this, this, new Identity<>(), new Identity<>(), new Identity<>(),
				new Identity<>());
	}

	@Override
	public LocalizeResult<S, S, S, S> localizeAtIdeal(Ideal<S> primeIdeal) {
		if (!primeIdeal.isPrime()) {
			throw new ArithmeticException("Not a prime ideal!");
		}
		return new LocalizeResult<>(this, primeIdeal, this, new Identity<>(), new Identity<>(),
				new ConstantMap<>(one()), new Identity<>());
	}

	@Override
	public Ideal<S> maximalIdeal() {
		return getZeroIdeal();
	}

	@Override
	public Field<S> reduction() {
		return this;
	}

	@Override
	public S reduce(S t) {
		return t;
	}

	@Override
	public S lift(S t) {
		return t;
	}

	@Override
	public boolean isInteger(S t) {
		return true;
	}

	@Override
	public S asInteger(S t) {
		return t;
	}
	
	@Override
	public S getDenominator(S t) {
		return t;
	}

	@Override
	public Field<S> reduction(Ideal<S> maximalIdeal) {
		return this;
	}

	@Override
	public S reduce(S t, Ideal<S> maximalIdeal) {
		return t;
	}

	@Override
	public S lift(S s, Ideal<S> maximalIdeal) {
		return s;
	}

	@Override
	public GlobalField<S, S, S> quotientField() {
		throw new ArithmeticException("implementation artifact");
	}

	@Override
	public DiscreteValuationRing<S, S> localize(Ideal<S> maximalIdeal) {
		return new LocalRingImplementation<>(new TrivialDiscreteValuationField<>(this), toString());
	}

	@Override
	public DiscreteValuationField<S, S> localizeAndQuotient(Ideal<S> maximalIdeal) {
		return new TrivialDiscreteValuationField<>(this);
	}

	@Override
	public boolean isZeroDivisor(S t) {
		return t.equals(zero());
	}

	@Override
	public boolean isEuclidean() {
		return true;
	}

	@Override
	public boolean isDivisible(S dividend, S divisor) {
		if (!divisor.equals(zero())) {
			return true;
		}
		return dividend.equals(zero());
	}

	@Override
	public PivotStrategy<S> preferredPivotStrategy() {
		return new TrivialPivotStrategy<>(this);
	}

	@Override
	public QuotientAndRemainderResult<S> quotientAndRemainder(S dividend, S divisor) {
		if (dividend.equals(zero()) && divisor.equals(zero())) {
			return new QuotientAndRemainderResult<>(zero(), zero());
		}
		return new QuotientAndRemainderResult<>(this.divide(dividend, divisor), zero());
	}

	@Override
	public BigInteger euclidMeasure(S t) {
		return BigInteger.ZERO;
	}

	public ExtendedEuclideanResult<S> extendedEuclidean(S t1, S t2) {
		if (t1.equals(this.zero()) && t2.equals(this.zero())) {
			return new ExtendedEuclideanResult<>(zero(), zero(), zero());
		}
		if (t1.equals(this.zero())) {
			return new ExtendedEuclideanResult<>(one(), zero(), inverse(t2));
		}
		if (t2.equals(this.zero()) || this.characteristic().equals(BigInteger.TWO)) {
			return new ExtendedEuclideanResult<>(one(), inverse(t1), zero());
		}
		return new ExtendedEuclideanResult<>(one(), this.inverse(this.multiply(2, t1)),
				this.inverse(this.multiply(2, t2)));
	}

	@Override
	public BezoutIdentityResult<S> bezoutIdentity(Ideal<S> t1, Ideal<S> t2) {
		if (!t1.contains(one()) && !t2.contains(one())) {
			throw new ArithmeticException("Ideals not coprime!");
		}
		if (!t1.contains(one())) {
			return new BezoutIdentityResult<>(Collections.singletonList(zero()), Collections.singletonList(one()));
		}
		if (!t2.contains(one()) || characteristic().equals(BigInteger.TWO)) {
			return new BezoutIdentityResult<>(Collections.singletonList(one()), Collections.singletonList(zero()));
		}
		S oneHalf = inverse(getInteger(2));
		return new BezoutIdentityResult<>(Collections.singletonList(oneHalf), Collections.singletonList(oneHalf));
	}

	@Override
	public Pair<S, S> bezoutIdentity(S t1, S t2) {
		if (t1.equals(zero()) && t2.equals(zero())) {
			throw new ArithmeticException("Elements not coprime!");
		}
		if (t1.equals(zero())) {
			return new Pair<>(zero(), inverse(t2));
		}
		if (t2.equals(zero()) || characteristic().equals(BigInteger.TWO)) {
			return new Pair<>(inverse(t1), zero());
		}
		return new Pair<>(inverse(multiply(2, t1)), inverse(multiply(2, t2)));
	}

	@Override
	public boolean coprime(Ideal<S> t1, Ideal<S> t2) {
		return t1.contains(one()) || t2.contains(one());
	}

	@Override
	public boolean coprime(S t1, S t2) {
		return !t1.equals(zero()) || !t2.equals(zero());
	}

	@Override
	public ChineseRemainderPreparation<S> prepareChineseRemainderTheorem(List<? extends Ideal<S>> ideals) {
		if (ideals.size() > 1) {
			throw new ArithmeticException("Not coprime proper ideals!");
		}
		if (ideals.size() == 1 && ideals.get(0).contains(one())) {
			throw new ArithmeticException("Not proper ideals!");
		}
		if (ideals.size() == 1) {
			return new ChineseRemainderPreparation<>(ideals, getZeroIdeal(), Collections.singletonList(one()));
		}
		return new ChineseRemainderPreparation<>(ideals, getZeroIdeal(), Collections.emptyList());
	}

	@Override
	public S chineseRemainderTheorem(List<S> elements, ChineseRemainderPreparation<S> preparation) {
		if (elements.size() != preparation.getMultipliers().size()) {
			throw new ArithmeticException("Mismatched multipliers!");
		}
		if (elements.size() == 0) {
			return zero();
		}
		return elements.get(0);
	}

	@Override
	public S chineseRemainderTheorem(List<S> elements, List<? extends Ideal<S>> ideals) {
		if (elements.size() != ideals.size()) {
			throw new ArithmeticException("Mismatched multipliers!");
		}
		if (elements.size() == 0) {
			return zero();
		}
		return elements.get(0);
	}

	@Override
	public Iterable<S> getUnits() {
		return getNonZeroElements();
	}
	
	/*@Override
	public FieldExtensionUnivariatePolynomialRing<T, S, Ext> getUnivariatePolynomialRing() {
		if (univariatePolynomialRing == null) {
			univariatePolynomialRing = new FieldExtensionUnivariatePolynomialRing<>(asExtensionType());
		}
		return univariatePolynomialRing;
	}*/

	@Override
	public FactorizationResult<Polynomial<S>, S> factorization(UnivariatePolynomial<S> t) {
		UnivariatePolynomialRing<S> ring = getUnivariatePolynomialRing();
		SortedMap<Polynomial<S>, Integer> result = new TreeMap<>();
		S unit = t.leadingCoefficient();
		t = ring.normalize(t);
		FactorizationResult<Polynomial<S>, S> squareFree = ring.squareFreeFactorization(t);
		for (Polynomial<S> squareFreeFactor : squareFree.primeFactors()) {
			for (Polynomial<S> factor : factorizeSquareFree(ring.toUnivariate(squareFreeFactor))) {
				result.put(factor, squareFree.multiplicity(squareFreeFactor));
			}
		}
		return new FactorizationResult<>(unit, result);
	}

	@SuppressWarnings("unchecked")
	private List<Polynomial<S>> factorizeSquareFree(UnivariatePolynomial<S> t) {
		UnivariatePolynomialRing<S> ring = getUnivariatePolynomialRing();
		if (t.degree() == 1) {
			return Collections.singletonList(ring.normalize(t));
		}
		Integers z = Integers.z();
		UnivariatePolynomialRing<T> rationalRing = baseField.getUnivariatePolynomialRing();
		GenericAlgebraicRingExtension<S> cr = new GenericAlgebraicRingExtension<>(t, this);
		Iterator<UnivariatePolynomial<IntE>> it = z.getUnivariatePolynomialRing().polynomials(degree() - 1);
		while (true) {
			S probeTerm = fromSmallDegreePolynomial(rationalRing.getEmbedding(it.next(), new MathMap<>() {
				@Override
				public T evaluate(IntE t) {
					return baseField.getInteger(t);
				}
			}));
			GenericAlgebraicExtensionElement<S> probe = cr.add(cr.alpha(),
					cr.getEmbedding(multiply(probeTerm, alpha())));
			List<Vector<T>> generators = new ArrayList<>();
			GenericAlgebraicExtensionElement<S> power = cr.one();
			for (int i = 0; i < t.degree() * degree(); i++) {
				generators.add(asRationalVector(power, t.degree()));
				power = cr.multiply(power, probe);
			}
			FiniteVectorSpace<T> asVectorSpace = new FiniteVectorSpace<>(baseField, t.degree() * degree());
			Matrix<T> m = Matrix.fromColumns(generators);
			MatrixAlgebra<T> algebra = asVectorSpace.matrixAlgebra();
			if (algebra.rank(m) != t.degree() * degree()) {
				continue;
			}
			UnivariatePolynomial<T> x = rationalRing.getPolynomial(
					algebra.solve(m, asRationalVector(cr.fromPolynomial(ring.getVar()), t.degree())).asList());
			UnivariatePolynomial<T> gamma = rationalRing
					.getPolynomial(algebra.solve(m, asRationalVector(cr.getEmbedding(alpha()), t.degree())).asList());
			Vector<T> lastPower = asRationalVector(power, t.degree());
			UnivariatePolynomial<T> minimalPolynomial = rationalRing
					.getPolynomial(algebra.solve(m, lastPower).asList());
			minimalPolynomial = rationalRing.toUnivariate(
					rationalRing.subtract(rationalRing.getVarPower(t.degree() * degree()), minimalPolynomial));
			FactorizationResult<Polynomial<T>, T> rationalFactors = baseField.factorization(minimalPolynomial);
			if (rationalFactors.primeFactors().size() == 1) {
				return Collections.singletonList(ring.normalize(t));
			}
			List<Polynomial<S>> factors = new ArrayList<>();
			for (Polynomial<T> factor : rationalFactors.primeFactors()) {
				Ext extension = makeExtension(rationalRing.toUnivariate(factor));
				if (!extension.minimalPolynomial().equals(factor)) {
					throw new ArithmeticException("Something messed with the minimalpolynomial!");
				}
				UnivariatePolynomialRing<S> extensionRing = extension.getUnivariatePolynomialRing();
				S embeddedX = extensionRing.evaluate(extensionRing.getEmbedding(x, extension.getEmbeddingMap()),
						extension.alpha());
				S embeddedGamma = extensionRing.evaluate(extensionRing.getEmbedding(gamma, extension.getEmbeddingMap()),
						extension.alpha());
				List<Vector<T>> gammaXBase = new ArrayList<>();
				S xPower = extension.one();
				for (int i = 0; i < extension.degree() / degree(); i++) {
					S gammaPower = extension.one();
					for (int j = 0; j < degree(); j++) {
						gammaXBase.add(extension.asVector(extension.multiply(xPower, gammaPower)));
						gammaPower = extension.multiply(gammaPower, embeddedGamma);
					}
					xPower = extension.multiply(xPower, embeddedX);
				}
				Matrix<T> gammaXMatrix = Matrix.fromColumns(gammaXBase);
				UnivariatePolynomial<S> extensionFactor = minimalPolynomial(embeddedX, extension.degree() / degree(),
						extension, this, new MathMap<>() {

							@Override
							public Vector<S> evaluate(S t) {
								Vector<T> overGammaX = extension.matrixAlgebra().solve(gammaXMatrix,
										extension.asVector(t));
								List<S> asList = new ArrayList<>();
								for (int i = 0; i < extension.degree() / degree(); i++) {
									asList.add(fromPolynomial(rationalRing.getPolynomial(
											overGammaX.asList().subList(i * degree(), (i + 1) * degree()))));
								}
								return new Vector<>(asList);
							}
						});
				factors.add(ring.normalize(extensionFactor));
			}
			return factors;
		}
	}

	private Vector<T> asRationalVector(GenericAlgebraicExtensionElement<S> e, int dimension) {
		UnivariatePolynomialRing<S> ring = getUnivariatePolynomialRing();
		List<S> asExtensionList = ring.asVector(e.asPolynomial(), dimension - 1).asList();
		List<T> asList = new ArrayList<>();
		for (S nfe : asExtensionList) {
			asList.addAll(asVector(nfe).asList());
		}
		return new Vector<>(asList);
	}

}
