package fields.polynomials;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.AbstractIdeal;
import fields.interfaces.Element;
import fields.interfaces.Ideal;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring;
import fields.interfaces.Ring.FactorizationResult;
import fields.interfaces.Ring.IdealResult;
import util.MiscAlgorithms;

public class PolynomialIdeal<T extends Element<T>> extends AbstractIdeal<Polynomial<T>> {
	private List<Polynomial<T>> basis;
	private List<List<Polynomial<T>>> relations;
	private PolynomialRing<T> polynomialRing;
	private CoordinateRing<T> coordinateRing;

	PolynomialIdeal(PolynomialRing<T> polynomialRing, List<Polynomial<T>> generators) {
		super(polynomialRing);
		this.polynomialRing = polynomialRing;
		this.basis = generators;
	}

	public CoordinateRing<T> divideOut() {
		if (coordinateRing == null) {
			coordinateRing = new CoordinateRing<>(polynomialRing, this);
		}
		return coordinateRing;
	}

	public int dimension() {
		return divideOut().krullDimension() - polynomialRing.getRing().krullDimension();
	}

	public int degree() {
		return divideOut().degree();
	}

	@Override
	public List<List<Polynomial<T>>> nonTrivialCombinations(List<Polynomial<T>> s) {
		Ring<T> ring = polynomialRing.getRing();
		IdealResult<Polynomial<T>, PolynomialIdeal<T>> ideal = polynomialRing.getIdealWithTransforms(s);
		if (ideal.getIdeal().generators().equals(s)) {
			return ideal.getIdeal().getModuleGeneratorRelations();
		}
		List<List<Polynomial<T>>> result = new ArrayList<>();
		for (int i = 1; i < s.size(); i++) {
			for (int j = 0; j < i; j++) {
				Monomial leadingFirst = s.get(j).leadingMonomial();
				T firstLeadingCoefficient = s.get(j).leadingCoefficient();
				Monomial leadingSecond = s.get(i).leadingMonomial();
				T secondLeadingCoefficient = s.get(i).leadingCoefficient();
				Monomial lcm = leadingFirst.lcm(leadingSecond);
				T coefficientLcm = ring.lcm(firstLeadingCoefficient, secondLeadingCoefficient);
				Polynomial<T> firstMultiplier = polynomialRing.getPolynomial(Collections.singletonMap(
						lcm.divide(leadingFirst), ring.divideChecked(coefficientLcm, firstLeadingCoefficient)));
				Polynomial<T> secondMultiplier = polynomialRing
						.getPolynomial(Collections.singletonMap(lcm.divide(leadingSecond),
								ring.negative(ring.divideChecked(coefficientLcm, secondLeadingCoefficient))));
				Polynomial<T> buchbergerPolynomial = polynomialRing.add(
						polynomialRing.multiply(firstMultiplier, s.get(j)),
						polynomialRing.multiply(secondMultiplier, s.get(i)));
				List<Polynomial<T>> generated = ideal.getIdeal().generate(buchbergerPolynomial);
				List<Polynomial<T>> relation = new ArrayList<>();
				for (int k = 0; k < s.size(); k++) {
					if (k == j) {
						relation.add(firstMultiplier);
					} else if (k == i) {
						relation.add(secondMultiplier);
					} else {
						relation.add(polynomialRing.zero());
					}
				}
				for (int k = 0; k < generated.size(); k++) {
					Polynomial<T> coefficient = polynomialRing.negative(generated.get(k));
					for (int l = 0; l < s.size(); l++) {
						Polynomial<T> originalCoefficient = polynomialRing.multiply(coefficient,
								ideal.getGeneratorExpressions().get(k).get(l));
						relation.set(l, polynomialRing.add(relation.get(l), originalCoefficient));
					}
				}
				result.add(relation);
			}
		}
		return result;
	}

	@Override
	public List<List<Polynomial<T>>> getModuleGeneratorRelations() {
		if (relations == null) {
			Ring<T> ring = polynomialRing.getRing();
			relations = new ArrayList<>();
			for (int i = 1; i < basis.size(); i++) {
				Monomial leadingFirst = basis.get(i - 1).leadingMonomial();
				T firstLeadingCoefficient = basis.get(i - 1).leadingCoefficient();
				Monomial leadingSecond = basis.get(i).leadingMonomial();
				T secondLeadingCoefficient = basis.get(i).leadingCoefficient();
				Monomial lcm = leadingFirst.lcm(leadingSecond);
				T coefficientLcm = ring.lcm(firstLeadingCoefficient, secondLeadingCoefficient);
				Polynomial<T> firstMultiplier = polynomialRing.getPolynomial(Collections.singletonMap(
						lcm.divide(leadingFirst), ring.divideChecked(coefficientLcm, firstLeadingCoefficient)));
				Polynomial<T> secondMultiplier = polynomialRing
						.getPolynomial(Collections.singletonMap(lcm.divide(leadingSecond),
								ring.negative(ring.divideChecked(coefficientLcm, secondLeadingCoefficient))));
				Polynomial<T> buchbergerPolynomial = polynomialRing.add(
						polynomialRing.multiply(firstMultiplier, basis.get(i - 1)),
						polynomialRing.multiply(secondMultiplier, basis.get(i)));
				List<Polynomial<T>> generated = generate(buchbergerPolynomial);
				List<Polynomial<T>> relation = new ArrayList<>();
				for (int j = 0; j < basis.size(); j++) {
					Polynomial<T> coefficient = polynomialRing.negative(generated.get(j));
					if (j == i - 1) {
						coefficient = polynomialRing.add(coefficient, firstMultiplier);
					} else if (j == i) {
						coefficient = polynomialRing.add(coefficient, secondMultiplier);
					}
					relation.add(coefficient);
				}
				relations.add(relation);
			}
		}
		return relations;
	}

	@Override
	public boolean contains(Polynomial<T> polynomial) {
		return polynomialRing.reduce(polynomial, this.basis).equals(polynomialRing.zero());
	}

	@Override
	public Polynomial<T> residue(Polynomial<T> f) {
		return polynomialRing.reduce(f, this.basis);
	}

	public PolynomialIdeal<T> saturate(Polynomial<T> by) {
		PolynomialRing<T> ring = polynomialRing.addVariableWithElimination(1);
		Polynomial<T> t = ring.getVar(1);
		Polynomial<T> f = ring.getEmbeddingWithElimination(by, 1);
		Polynomial<T> inv = ring.subtract(ring.one(), ring.multiply(f, t));
		List<Polynomial<T>> generators = new ArrayList<Polynomial<T>>();
		generators.add(inv);
		for (Polynomial<T> p : this.basis) {
			generators.add(ring.getEmbeddingWithElimination(p, 1));
		}
		Ideal<Polynomial<T>> saturationIdeal = ring.getIdeal(generators);
		List<Polynomial<T>> saturationGenerators = new ArrayList<Polynomial<T>>();
		for (Polynomial<T> b : saturationIdeal.generators()) {
			if (b.leadingMonomial().exponents()[0] == 0)
				saturationGenerators.add(polynomialRing.getEmbeddingWithElimination(b, -1));
		}
		return new PolynomialIdeal<T>(polynomialRing, saturationGenerators);

	}

	public String toString() {
		return this.basis.toString();
	}

	public boolean isFinite() {
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() {
		throw new InfinityException();
	}

	@Override
	public List<Polynomial<T>> generators() {
		return Collections.unmodifiableList(this.basis);
	}

	@Override
	public List<Polynomial<T>> generate(Polynomial<T> t) {
		return polynomialRing.generalQuotientAndRemainder(t, basis).getQuotients();
	}

	public Ideal<T> intersectToRing() {
		List<T> ringGenerators = new ArrayList<>();
		for (Polynomial<T> generator : generators()) {
			if (generator.degree() == 0) {
				ringGenerators.add(generator.leadingCoefficient());
			}
		}
		return polynomialRing.getRing().getIdeal(ringGenerators);
	}

	public boolean isHomogenous() {
		for (Polynomial<T> generator : generators()) {
			if (!polynomialRing.isHomogeneous(generator)) {
				return false;
			}
		}
		return true;
	}

	@Override
	public boolean isPrimary() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean isPrime() {
		if (!intersectToRing().isPrime()) {
			return false;
		}
		if (polynomialRing.numberOfVariables() == 0) {
			return true;
		}
		if (polynomialRing.numberOfVariables() == 1) {
			for (Polynomial<T> generator : generators()) {
				if (generator.degree() > 0 && !polynomialRing.isIrreducible(generator)) {
					return false;
				}
			}
			return true;
		}
		if (polynomialRing.getComparator() != Monomial.REVLEX) {
			PolynomialRing<T> revlexRing = AbstractPolynomialRing.getPolynomialRing(polynomialRing.getRing(),
					polynomialRing.numberOfVariables(), Monomial.REVLEX);
			return revlexRing.getEmbedding(this).isPrime();
		}
		List<Polynomial<T>> intersect = new ArrayList<>();
		for (Polynomial<T> polynomial : generators()) {
			if (polynomial.degree(polynomialRing.numberOfVariables()) == 0) {
				intersect.add(polynomialRing.eliminateVariable().getEmbedding(polynomial));
			}
		}
		PolynomialIdeal<T> intersectIdeal = polynomialRing.eliminateVariable().getIdeal(intersect);
		if (!intersectIdeal.isPrime()) {
			return false;
		}
		return intersectToRing().isPrime();
	}

	@Override
	public boolean isMaximal() {
		return isPrime() && divideOut().krullDimension() == polynomialRing.getRing().krullDimension()
				&& intersectToRing().isMaximal();
	}

	@Override
	public boolean isRadical() {
		for (Polynomial<T> polynomial : generators()) {
			if (!polynomialRing.isSquareFree(polynomial)) {
				return false;
			}
		}
		return intersectToRing().isRadical();
	}

	public List<PolynomialIdeal<T>> minimalPrimeIdealsOver() {
		if (!polynomialRing.getRing().isIntegral()) {
			throw new UnsupportedOperationException("Non integral base ring!");
		}
		if (!isRadical()) {
			return polynomialRing.radical(this).minimalPrimeIdealsOver();
		}
		List<List<Polynomial<T>>> factors = new ArrayList<>();
		for (Polynomial<T> generator : generators()) {
			if (generator.equals(polynomialRing.zero())) {
				continue;
			}
			FactorizationResult<Polynomial<T>, Polynomial<T>> factorization = polynomialRing
					.uniqueFactorization(generator);
			List<Polynomial<T>> factorsOfGenerator = new ArrayList<>();
			factorsOfGenerator.addAll(factorization.primeFactors());
			factors.add(factorsOfGenerator);
		}
		List<List<Polynomial<T>>> crossProduct = MiscAlgorithms.crossProduct(factors);
		List<PolynomialIdeal<T>> result = new ArrayList<>();
		for (List<Polynomial<T>> crossed : crossProduct) {
			result.add(polynomialRing.getIdeal(crossed));
		}
		return result;
	}

	public boolean isMaximalOverAlgebraicClosure() {
		if (polynomialRing.numberOfVariables() != polynomialRing.krullDimension()) {
			throw new ArithmeticException("Only defined over field!");
		}
		if (polynomialRing.numberOfVariables() != generators().size()) {
			return false;
		}
		for (Polynomial<T> generator : generators()) {
			if (generator.degree() != 1) {
				return false;
			}
		}
		return true;
	}

}
