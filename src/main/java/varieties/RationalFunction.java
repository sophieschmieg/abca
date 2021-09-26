package varieties;

import java.util.ArrayList;
import java.util.List;

import fields.helper.CoordinateRing;
import fields.interfaces.Element;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
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

	public RationalFunction(ProjectiveScheme<T> domain, ProjectiveLine<T> range,
			ProjectiveMorphism<T> morphism) {
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
		this.dehomogenizedNumerator = affinePolynomialRing.getEmbedding(
				domain.asGenericProjectiveScheme().homogenousPolynomialRing().dehomogenize(numerator, affineCoverIndex + 1));
		this.dehomogenizedDenominator = affinePolynomialRing.getEmbedding(domain.asGenericProjectiveScheme()
				.homogenousPolynomialRing().dehomogenize(denominator, affineCoverIndex + 1));
	}
	
	@Override
	public ProjectiveScheme<T> getDomain() {
		return domain;
	}

	@Override
	public ProjectiveLine<T> getRange() {
		return range;
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

	public Polynomial<T> getNumerator() {
		return numerator;
	}

	public Polynomial<T> getDenominator() {
		return denominator;
	}
}
