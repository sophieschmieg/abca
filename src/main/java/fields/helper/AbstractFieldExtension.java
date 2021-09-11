package fields.helper;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.GenericAlgebraicRingExtension.GenericAlgebraicExtensionElement;
import fields.integers.Rationals.Fraction;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.BilinearMap;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.FieldExtension;
import fields.interfaces.Ideal;
import fields.interfaces.MathMap;
import fields.interfaces.Module;
import fields.interfaces.Polynomial;
import fields.interfaces.Ring;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.vectors.FiniteVectorSpace;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;
import fields.vectors.pivot.PivotStrategy;
import fields.vectors.pivot.TrivialPivotStrategy;

public abstract class AbstractFieldExtension<T extends Element<T>, S extends AlgebraicExtensionElement<T, S>, Ext extends FieldExtension<T, S, Ext>>
		extends AbstractAlgebraicRingExtension<T, S, Ext> implements FieldExtension<T, S, Ext> {
	private Field<T> baseField;
	private UnivariatePolynomialRing<T> polynomials;
	private UnivariatePolynomial<T> minimalPolynomial;
	private int degree;
	private FiniteVectorSpace<T> asVectorSpace;
	private Polynomial<T> genericNorm;
	
	public AbstractFieldExtension(UnivariatePolynomial<T> minimalPolynomial, Field<T> baseField) {
		super(minimalPolynomial, baseField);
		this.baseField = baseField;
		this.polynomials = baseField.getUnivariatePolynomialRing();
		this.minimalPolynomial = minimalPolynomial;
		this.degree = this.minimalPolynomial.degree();
		this.asVectorSpace = new FiniteVectorSpace<>(this.baseField, degree);
	}

	public AbstractFieldExtension(Field<T> baseField) {
		this(baseField.getUnivariatePolynomialRing().getVar(), baseField);
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
					asVector.addAll(AbstractFieldExtension.this.asVector(t.asPolynomial().univariateCoefficient(i)).asList());
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
		FieldEmbedding<T, S, Ext> extension = getEmbeddedExtension(minimalPolynomial);
		return new Extension<S, T, S, Ext>(extension.getField(), this, extension.getEmbeddingMap(), new MathMap<>() {
			@Override
			public Vector<S> evaluate(S t) {
				return extension.asVector(t);
			}
		});
	}

	@Override
	public boolean isIrreducible(UnivariatePolynomial<S> t) {
		return factorization(t).isIrreducible();
	}

	@Override
	public int krullDimension() {
		return 0;
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
		throw new ArithmeticException("Ideal not proper and non zero (it's a field!)");
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

	@Override
	public Iterable<S> getUnits() {
		return getNonZeroElements();
	}

}
