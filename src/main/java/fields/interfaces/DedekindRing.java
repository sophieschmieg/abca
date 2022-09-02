package fields.interfaces;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import fields.interfaces.DiscreteValuationRing.OkutsuType;
import fields.interfaces.DiscreteValuationRing.TheMontesResult;
import fields.interfaces.Field.Extension;
import fields.local.Value;
import util.FunctionMathMap;

public interface DedekindRing<T extends Element<T>, U extends Element<U>, S extends Element<S>> extends Ring<T> {
	Value valuation(T t, Ideal<T> maximalIdeal);

	FieldOfFractionsResult<T, U> fieldOfFractions();

	boolean isInteger(U t);

	@Override
	default boolean isIntegral() {
		return true;
	}

	@Override
	default boolean isReduced() {
		return true;
	}

	@Override
	default boolean isIrreducible() {
		return true;
	}

	@Override
	default boolean isDedekindDomain() {
		return true;
	}

	default DedekindRing<T, U, S> asDedekindRing() {
		return this;
	}

	T asInteger(U t);

	T getDenominator(U t);

	default T getNumerator(U t) {
		return asInteger(quotientField().multiply(quotientField().getInteger(getDenominator(t)), t));
	}

	default T uniformizer(Ideal<T> maximalIdeal) {
		return asInteger(localize(maximalIdeal).uniformizer());
	}

	public static class TwoGeneratorIdeal<T extends Element<T>, U extends Element<U>, S extends Element<S>, R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>> {
		private UnivariatePolynomial<T> minimalPolynomial;
		private UnivariatePolynomial<U> uniformizer;
		private T intGenerator;
		private OkutsuType<U, S, R, RE, RFE> type;
		private int index;
		private List<OkutsuType<U, S, R, RE, RFE>> types;
		private List<List<UnivariatePolynomial<U>>> integralBases;
		private DiscreteValuationRing<U, S> localRing;

		private TwoGeneratorIdeal(UnivariatePolynomial<T> minimalPolynomial, UnivariatePolynomial<U> uniformizer,
				T intGenerator, OkutsuType<U, S, R, RE, RFE> type, int index, List<OkutsuType<U, S, R, RE, RFE>> types,
				List<List<UnivariatePolynomial<U>>> integralBases, DiscreteValuationRing<U, S> localRing) {
			this.minimalPolynomial = minimalPolynomial;
			this.uniformizer = uniformizer;
			this.intGenerator = intGenerator;
			this.type = type;
			this.index = index;
			this.types = types;
			this.integralBases = integralBases;
			this.localRing = localRing;
		}

		public UnivariatePolynomial<T> getMinimalPolynomial() {
			return minimalPolynomial;
		}

		public UnivariatePolynomial<U> getUniformizer() {
			return uniformizer;
		}

		public T getIntGenerator() {
			return intGenerator;
		}

		public OkutsuType<U, S, R, RE, RFE> getType() {
			return type;
		}

		public int getIndex() {
			return index;
		}

		public List<OkutsuType<U, S, R, RE, RFE>> getTypes() {
			return types;
		}

		public List<List<UnivariatePolynomial<U>>> getIntegralBases() {
			return integralBases;
		}

		public DiscreteValuationRing<U, S> getLocalRing() {
			return localRing;
		}
	}

	default <R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>> List<TwoGeneratorIdeal<T, U, S, R, RE, RFE>> idealsOverGenerators(
			UnivariatePolynomial<T> minimalPolynomial, T prime, Extension<S, R, RE, RFE> trivialExtension) {
		if (!isPrime(prime)) {
			throw new ArithmeticException("Not a maximal principal ideal!");
		}
		DiscreteValuationRing<U, S> local = localize(getIdeal(Collections.singletonList(prime)));
		GlobalField<U, T, S> field = quotientField();
		UnivariatePolynomial<U> localMinimalPolynomial = field.getUnivariatePolynomialRing()
				.getEmbedding(minimalPolynomial, new FunctionMathMap<>((T t) -> field.getInteger(t)));
		TheMontesResult<U, S, R, RE, RFE> theMontes = local.theMontesAlgorithm(localMinimalPolynomial,
				trivialExtension);
		List<TwoGeneratorIdeal<T, U, S, R, RE, RFE>> ideals = new ArrayList<>();
		List<List<UnivariatePolynomial<U>>> integralBasis = local.integralBasis(localMinimalPolynomial, theMontes,
				true);
		for (int index = 0; index < theMontes.getTypes().size(); index++) {
			OkutsuType<U, S, R, RE, RFE> type = theMontes.getTypes().get(index);
			UnivariatePolynomialRing<U> polynomials = field.getUnivariatePolynomialRing();
			UnivariatePolynomial<U> uniformizer = null;
			UnivariatePolynomial<U> basisElement = null;
			if (type.ramificationIndex() == 1) {
				uniformizer = polynomials.getEmbedding(local.uniformizer());
				basisElement = integralBasis.get(index).get(0);
			} else {
				boolean foundUniformizer = false;
				boolean foundZero = false;
				for (int i = 0; i < integralBasis.get(index).size(); i++) {
					Value value = type.valuation(integralBasis.get(index).get(i));
					if (!foundUniformizer && value.equals(Value.ONE)) {
						foundUniformizer = true;
						uniformizer = integralBasis.get(index).get(i);
					} else if (!foundZero && value.equals(Value.ZERO)) {
						foundZero = true;
						basisElement = integralBasis.get(index).get(i);
					}
					if (foundZero && foundUniformizer) {
						break;
					}
				}
				if (!foundUniformizer) {
					System.err.println("Did not find uniformizer, but kinda expected to find one!");
					uniformizer = type.lift(type.reduction().extension().one(), 1);
				}
				if (!foundZero) {
					throw new ArithmeticException("Did not find zero valuation element in integral basis!");
				}
			}
			UnivariatePolynomial<U> generator = polynomials.multiply(uniformizer, basisElement);
			for (int i = 0; i < theMontes.getTypes().size(); i++) {
				if (i == index) {
					continue;
				}
				OkutsuType<U, S, R, RE, RFE> otherType = theMontes.getTypes().get(i);
				if (!otherType.valuation(generator).equals(Value.ZERO)) {
					for (int j = 0; j < integralBasis.get(i).size(); j++) {
						if (otherType.valuation(integralBasis.get(i).get(j)).equals(Value.ZERO)) {
							generator = polynomials.add(generator, integralBasis.get(i).get(j));
							break;
						}
					}
				}
			}
			for (int i = 0; i < theMontes.getTypes().size(); i++) {
				if (i == index) {
					if (!theMontes.getTypes().get(i).valuation(generator).equals(Value.ONE)) {
						throw new ArithmeticException("Not a uniformizer of the prime ideal!");
					}
				} else {
					if (!theMontes.getTypes().get(i).valuation(generator).equals(Value.ZERO)) {
						throw new ArithmeticException("Not a unit mod the other prime ideals!");
					}
				}
			}
			ideals.add(new TwoGeneratorIdeal<T, U, S, R, RE, RFE>(minimalPolynomial, generator,
					asInteger(local.uniformizer()), type, index, theMontes.getTypes(), integralBasis, local));
		}
		return ideals;
	}

	default List<UnivariatePolynomial<U>> extensionIntegralBasis(UnivariatePolynomial<T> minimalPolynomial) {
		GlobalField<U, T, S> field = quotientField();
		UnivariatePolynomialRing<U> polynomials = field.getUnivariatePolynomialRing();
		if (minimalPolynomial.degree() == 1) {
			return Collections.singletonList(polynomials.one());
		}
		UnivariatePolynomialRing<T> integerPolynomials = getUnivariatePolynomialRing();
		List<T> primes = new ArrayList<>();
		primes.addAll(uniqueFactorization(integerPolynomials.discriminant(minimalPolynomial)).primeFactors());
		Map<T, List<UnivariatePolynomial<U>>> integralBasisPerPrime = new TreeMap<>();
		Map<T, List<Integer>> valuesPerPrime = new TreeMap<>();
		UnivariatePolynomial<U> embeddedMinimalPolynomial = polynomials.getEmbedding(minimalPolynomial,
				new FunctionMathMap<>((T t) -> field.getInteger(t)));
		for (T prime : primes) {
			DiscreteValuationRing<U, S> localized = localize(getIdeal(Collections.singletonList(prime)));
			List<UnivariatePolynomial<U>> integral = localized.triagonalizeIntegralBasis(embeddedMinimalPolynomial,
					localized.integralBasis(embeddedMinimalPolynomial,
							localized.theMontesAlgorithm(embeddedMinimalPolynomial), false));
			List<Integer> valueList = new ArrayList<>();
			for (int i = 0; i < minimalPolynomial.degree(); i++) {
				valueList.add(-localized.localField().valuation(integral.get(i).leadingCoefficient()).value());
			}
			integralBasisPerPrime.put(prime, integral);
			valuesPerPrime.put(prime, valueList);
		}
		List<UnivariatePolynomial<U>> integralBasis = new ArrayList<>();
		for (int i = 0; i < minimalPolynomial.degree(); i++) {
			UnivariatePolynomial<U> b = polynomials.zero();
			List<T> moduli = new ArrayList<>();
			T divisor = one();
			for (T prime : primes) {
				T primePower = power(prime, valuesPerPrime.get(prime).get(i));
				divisor = multiply(primePower, divisor);
				moduli.add(multiply(primePower, prime));
			}
			ChineseRemainderPreparation<T> preparation = prepareChineseRemainderTheoremModuli(moduli);
			for (int j = 0; j < moduli.size(); j++) {
				T prime = primes.get(j);
				T primePower = power(prime, valuesPerPrime.get(prime).get(i));
					b = polynomials.add(b,
						polynomials.multiply(field.getInteger(multiply(preparation.getMultipliers().get(j), primePower)),
								integralBasisPerPrime.get(prime).get(i)));
			}
			b = polynomials.getEmbedding(b,
					new FunctionMathMap<>((U t) -> field.getInteger(preparation.getProduct().residue(asInteger(t)))));
			b = polynomials.divideScalar(b, field.getInteger(divisor));
			if (b.degree() != i) {
				throw new ArithmeticException("expected integral basis to be triagonal!");
			}
			integralBasis.add(b);
		}
		return integralBasis;
	}

	Field<S> reduction(Ideal<T> maximalIdeal);

	S reduce(T t, Ideal<T> maximalIdeal);

	T lift(S s, Ideal<T> maximalIdeal);

	DiscreteValuationRing<U, S> localize(Ideal<T> maximalIdeal);

	DiscreteValuationField<U, S> localizeAndQuotient(Ideal<T> maximalIdeal);

	GlobalField<U, T, S> quotientField();
}
