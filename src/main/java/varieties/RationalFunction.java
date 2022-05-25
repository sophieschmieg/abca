package varieties;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import fields.interfaces.Element;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.CoordinateRing;
import varieties.affine.AffinePoint;
import varieties.affine.AffineScheme;
import varieties.curves.ProjectiveLine;
import varieties.projective.ProjectiveMorphism;
import varieties.projective.ProjectivePoint;
import varieties.projective.ProjectiveScheme;

public class RationalFunction<T extends Element<T>>
		implements Morphism<T, ProjectivePoint<T>, ProjectivePoint<T>>, Element<RationalFunction<T>> {
	private ProjectiveScheme<T> domain;
	private ProjectiveLine<T> range;
	private ProjectiveMorphism<T> asProjectiveMorphism;
	private Polynomial<T> numerator;
	private Polynomial<T> denominator;
	private AffineScheme<T> affineSlice;
	private CoordinateRing<T> affineCoordinateRing;
	private PolynomialRing<T> affinePolynomialRing;
	private int affineCoverIndex;
	private Polynomial<T> dehomogenizedNumerator;
	private Polynomial<T> dehomogenizedDenominator;

	public static <T extends Element<T>> RationalFunction<T> fromFraction(ProjectiveScheme<T> domain,
			ProjectiveLine<T> range, Polynomial<T> numerator, Polynomial<T> denominator) {
		List<Polynomial<T>> asList = new ArrayList<>();
		asList.add(numerator);
		asList.add(denominator);
		ProjectiveMorphism<T> asProjectiveMorphism = new ProjectiveMorphism<>(domain.asGenericProjectiveScheme(),
				range.asGenericProjectiveScheme(), asList);
		return new RationalFunction<>(domain, range, asProjectiveMorphism);
	}

	public RationalFunction(ProjectiveScheme<T> domain, ProjectiveLine<T> range, ProjectiveMorphism<T> morphism) {
		if (morphism.getRange().dimension() != 1) {
			throw new ArithmeticException("Not a meromorphic function!");
		}
		if (morphism.getRange().projectiveEmbeddingDimension() != 1) {
			throw new ArithmeticException("Not a meromorphic function!");
		}
		this.domain = domain;
		this.asProjectiveMorphism = morphism;
		this.range = range;
		this.numerator = asProjectiveMorphism.asPolynomials().get(0);
		this.denominator = asProjectiveMorphism.asPolynomials().get(1);
		this.affineCoverIndex = domain.asGenericProjectiveScheme().homogenousPolynomialRing().numberOfVariables() - 1;
		this.affineSlice = domain.getAffineCover().getCover().get(affineCoverIndex);
		this.affineCoordinateRing = affineSlice.getCoordinateRing();
		this.affinePolynomialRing = affineCoordinateRing.getPolynomialRing();
		this.dehomogenizedNumerator = affinePolynomialRing.getEmbedding(domain.asGenericProjectiveScheme()
				.homogenousPolynomialRing().dehomogenize(numerator, affineCoverIndex + 1));
		this.dehomogenizedDenominator = affinePolynomialRing.getEmbedding(domain.asGenericProjectiveScheme()
				.homogenousPolynomialRing().dehomogenize(denominator, affineCoverIndex + 1));
	}
	
	@Override
	public String toString() {
		return numerator + "/" + denominator;
	}

	@Override
	public ProjectiveScheme<T> getDomain() {
		return domain;
	}

	@Override
	public ProjectiveLine<T> getRange() {
		return range;
	}

	private List<ProjectivePoint<T>> findRoots(Polynomial<T> asHomogenousPolynomial,
			Polynomial<T> asDehomogenizedPolynomial, ProjectivePoint<T> target) {
		List<ProjectivePoint<T>> candidates = new ArrayList<>();
		List<AffinePoint<T>> affineSolutions = affinePolynomialRing.solve(affineCoordinateRing
				.getIdeal(Collections.singletonList(affineCoordinateRing.getEmbedding(asDehomogenizedPolynomial)))
				.asPolynomialIdeal());
		for (AffinePoint<T> point : affineSolutions) {
			List<T> projective = new ArrayList<>();
			for (int i = 0; i < affinePolynomialRing.numberOfVariables(); i++) {
				if (i == affineCoverIndex) {
					projective.add(domain.getField().one());
				}
				projective.add(point.getCoord(i + 1));
			}
			if (affineCoverIndex == affinePolynomialRing.numberOfVariables()) {
				projective.add(domain.getField().one());
			}
			candidates.add(new ProjectivePoint<>(domain.getField(), projective));
		}
		List<T> eval = new ArrayList<>();
		for (int i = 0; i < affinePolynomialRing.numberOfVariables(); i++) {
			if (i == affineCoverIndex) {
				eval.add(domain.getField().zero());
			}
			eval.add(null);
		}
		if (affineCoverIndex == affinePolynomialRing.numberOfVariables()) {
			eval.add(domain.getField().zero());
		}
		PolynomialRing<T> homogenous = domain.asGenericProjectiveScheme().homogenousPolynomialRing();
		List<Polynomial<T>> infinityPoints = new ArrayList<>();
		infinityPoints.add(affinePolynomialRing.getEmbedding(homogenous.partiallyEvaluate(asHomogenousPolynomial, eval)));
		for (Polynomial<T> generator : domain.asGenericProjectiveScheme().generators()) {
			infinityPoints.add(affinePolynomialRing.getEmbedding(homogenous.partiallyEvaluate(generator, eval)));
		}
		List<AffinePoint<T>> infinitySolutions = affinePolynomialRing.solve(infinityPoints);
		for (AffinePoint<T> point : infinitySolutions) {
			List<T> projective = new ArrayList<>();
			boolean nonZero = false;
			for (int i = 0; i < affinePolynomialRing.numberOfVariables(); i++) {
				if (i == affineCoverIndex) {
					projective.add(domain.getField().zero());
				}
				projective.add(point.getCoord(i + 1));
				if (!nonZero && !point.getCoord(i + 1).equals(domain.getField().zero())) {
					nonZero = true;
				}
			}
			if (nonZero) {
				candidates.add(new ProjectivePoint<>(domain.getField(), projective));
			}
		}
		List<ProjectivePoint<T>> result = new ArrayList<>();
		for (ProjectivePoint<T> candidate : candidates) {
			if (evaluate(candidate).equals(target)) {
				result.add(candidate);
			}
		}
		return result;
	}

	public List<ProjectivePoint<T>> getZeroes() {
		return findRoots(numerator, dehomogenizedNumerator,
				new ProjectivePoint<>(domain.getField(), domain.getField().zero(), domain.getField().one()));
	}

	public List<ProjectivePoint<T>> getPoles() {
		return findRoots(denominator, dehomogenizedDenominator,
				new ProjectivePoint<>(domain.getField(), domain.getField().one(), domain.getField().zero()));
	}

	@Override
	public ProjectivePoint<T> evaluate(ProjectivePoint<T> t) {
		return asProjectiveMorphism.evaluate(t);
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof RationalFunction<?>)) {
			return false;
		}
		@SuppressWarnings("unchecked")
		RationalFunction<T> o = (RationalFunction<T>) obj;
		Polynomial<T> lhs = affinePolynomialRing.multiply(dehomogenizedNumerator, o.dehomogenizedDenominator);
		Polynomial<T> rhs = affinePolynomialRing.multiply(o.dehomogenizedNumerator, dehomogenizedDenominator);
		return affineCoordinateRing.getEmbedding(lhs).equals(affineCoordinateRing.getEmbedding(rhs));
	}

	@Override
	public int compareTo(RationalFunction<T> o) {
		Polynomial<T> gcd = affinePolynomialRing.gcd(dehomogenizedNumerator, dehomogenizedDenominator);
		Polynomial<T> reducedNumerator = affinePolynomialRing.divideChecked(dehomogenizedNumerator, gcd);
		Polynomial<T> reducedDenominator = affinePolynomialRing.divideChecked(dehomogenizedDenominator, gcd);
		if (reducedDenominator.compareTo(affinePolynomialRing.zero()) < 0) {
			reducedNumerator = affinePolynomialRing.multiply(-1, reducedNumerator);
			reducedDenominator = affinePolynomialRing.multiply(-1, reducedDenominator);
		}
		Polynomial<T> otherGcd = affinePolynomialRing.gcd(o.dehomogenizedNumerator, o.dehomogenizedDenominator);
		Polynomial<T> otherReducedNumerator = affinePolynomialRing.divideChecked(o.dehomogenizedNumerator, otherGcd);
		Polynomial<T> otherReducedDenominator = affinePolynomialRing.divideChecked(o.dehomogenizedDenominator,
				otherGcd);
		if (otherReducedDenominator.compareTo(affinePolynomialRing.zero()) < 0) {
			otherReducedNumerator = affinePolynomialRing.multiply(-1, otherReducedNumerator);
			otherReducedDenominator = affinePolynomialRing.multiply(-1, otherReducedDenominator);
		}
		Polynomial<T> lhs = affinePolynomialRing.multiply(reducedNumerator, otherReducedDenominator);
		Polynomial<T> rhs = affinePolynomialRing.multiply(otherReducedNumerator, reducedDenominator);
		return affineCoordinateRing.getEmbedding(lhs).compareTo(affineCoordinateRing.getEmbedding(rhs));
	}

	@Override
	public RestrictionResult<T> restrict(ProjectivePoint<T> preimage) {
		return asProjectiveMorphism.restrict(preimage);
	}

	public ProjectiveMorphism<T> asProjectiveMorphism() {
		return asProjectiveMorphism;
	}

	public Polynomial<T> getNumerator() {
		return numerator;
	}

	public Polynomial<T> getDenominator() {
		return denominator;
	}
}
