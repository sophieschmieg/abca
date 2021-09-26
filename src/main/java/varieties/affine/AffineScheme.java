package varieties.affine;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.CoordinateRing;
import fields.helper.CoordinateRing.CoordinateIdeal;
import fields.helper.CoordinateRing.CoordinateRingElement;
import fields.helper.LocalizedCoordinateRing;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.PolynomialIdeal;
import varieties.AbstractScheme;
import varieties.Morphism;
import varieties.SpectrumOfField;
import varieties.SpectrumOfField.SingletonPoint;

public class AffineScheme<T extends Element<T>> extends AbstractScheme<T, AffinePoint<T>> {
	private CoordinateRing<T> coordinateRing;
	private Field<T> field;
	private AffineCover<T> cover;

	public AffineScheme(Field<T> field, CoordinateRing<T> coordinateRing) {
		if (!coordinateRing.getPolynomialRing().getRing().equals(field)) {
			throw new ArithmeticException("Incompatible parameters");
		}
		this.coordinateRing = coordinateRing;
		this.field = field;
		this.cover = new AffineCover<>(Collections.singletonList(this),
				Collections.singletonList(Collections.singletonList(identityMorphism())));
	}

	public static class IntersectionResult<T extends Element<T>> {
		private AffineScheme<T> intersection;
		private AffineMorphism<T> firstEmbedding;
		private AffineMorphism<T> secondEmbedding;

		public AffineScheme<T> getIntersection() {
			return intersection;
		}

		public AffineMorphism<T> getFirstEmbedding() {
			return firstEmbedding;
		}

		public AffineMorphism<T> getSecondEmbedding() {
			return secondEmbedding;
		}
	}

	public static <T extends Element<T>> IntersectionResult<T> intersect(AffineScheme<T> t1, AffineScheme<T> t2) {
		PolynomialRing<T> polynomialRing = t1.coordinateRing.getPolynomialRing();
		if (!polynomialRing.equals(t2.coordinateRing.getPolynomialRing())) {
			throw new ArithmeticException("different polynomial rings!");
		}
		PolynomialIdeal<T> ideal = polynomialRing.add(t1.coordinateRing.getIdeal(), t2.coordinateRing.getIdeal());
		List<Polynomial<T>> identity = new ArrayList<>();
		for (int i = 0; i < polynomialRing.numberOfVariables(); i++) {
			identity.add(polynomialRing.getVar(i + 1));
		}
		IntersectionResult<T> result = new IntersectionResult<>();
		result.intersection = new AffineScheme<>(t1.getField(), new CoordinateRing<>(polynomialRing, ideal));
		result.firstEmbedding = AffineMorphism.fromPolynomials(result.intersection, t1, identity);
		result.secondEmbedding = AffineMorphism.fromPolynomials(result.intersection, t2, identity);
		return result;
	}

	public static <T extends Element<T>> AffineMorphism<T> restrictAwayFrom(AffineScheme<T> variety,
			List<CoordinateRingElement<T>> constraints) {
		int numVars = variety.coordinateRing.getPolynomialRing().numberOfVariables();
		PolynomialRing<T> polynomialRing = AbstractPolynomialRing.getPolynomialRing(variety.field, numVars + constraints.size(), variety.coordinateRing.getPolynomialRing().getComparator());
		List<Polynomial<T>> idealGenerators = new ArrayList<>();
		for(Polynomial<T> generator : variety.coordinateRing.getIdeal().generators()) {
			idealGenerators.add(polynomialRing.getEmbedding(generator));
		}
		for (int i = 0; i < constraints.size();i++) {
			Polynomial<T> t = polynomialRing.getEmbedding(constraints.get(i).getElement());
			idealGenerators.add(polynomialRing.subtract(polynomialRing.multiply(t, polynomialRing.getVar(i+1+numVars)), polynomialRing.one()));
		}
		PolynomialIdeal<T> ideal = polynomialRing.getIdeal(idealGenerators);
		List<Polynomial<T>> identity = new ArrayList<>();
		for (int i = 0; i < numVars; i++) {
			identity.add(polynomialRing.getVar(i + 1));
		}
		AffineScheme<T> restricted = new AffineScheme<>(variety.getField(),
				new CoordinateRing<T>(polynomialRing, ideal));
		return AffineMorphism.fromPolynomials(restricted, variety, identity);
	}

	public static class UnionResult<T extends Element<T>> {
		private AffineScheme<T> union;
		private AffineMorphism<T> firstEmbedding;
		private AffineMorphism<T> secondEmbedding;

		public AffineScheme<T> getUnion() {
			return union;
		}

		public AffineMorphism<T> getFirstEmbedding() {
			return firstEmbedding;
		}

		public AffineMorphism<T> getSecondEmbedding() {
			return secondEmbedding;
		}
	}

	public static <T extends Element<T>> UnionResult<T> union(AffineScheme<T> t1, AffineScheme<T> t2) {
		PolynomialRing<T> polynomialRing = t1.coordinateRing.getPolynomialRing();
		if (!polynomialRing.equals(t2.coordinateRing.getPolynomialRing())) {
			throw new ArithmeticException("different polynomial rings!");
		}
		PolynomialIdeal<T> ideal = polynomialRing.intersect(t1.coordinateRing.getIdeal(), t2.coordinateRing.getIdeal());
		List<Polynomial<T>> identity = new ArrayList<>();
		for (int i = 0; i < polynomialRing.numberOfVariables(); i++) {
			identity.add(polynomialRing.getVar(i + 1));
		}
		UnionResult<T> result = new UnionResult<>();
		result.union = new AffineScheme<>(t1.getField(), new CoordinateRing<>(polynomialRing, ideal));
		result.firstEmbedding = AffineMorphism.fromPolynomials(t1, result.union, identity);
		result.secondEmbedding = AffineMorphism.fromPolynomials(t2, result.union, identity);
		return result;
	}

	@Override
	public String toString() {
		return "Spec " + coordinateRing;
	}

	public CoordinateRing<T> getCoordinateRing() {
		return coordinateRing;
	}

	public LocalizedCoordinateRing<T> localizedCoordinateRing(AffinePoint<T> point) {
		CoordinateIdeal<T> ideal = coordinateRing.getIdeal(point.asIdeal(coordinateRing.getPolynomialRing()));
		return new LocalizedCoordinateRing<>(field, coordinateRing, ideal);
	}

	public boolean hasRationalPoint(AffinePoint<T> point) {
		return point.asIdeal(coordinateRing.getPolynomialRing()).contains(coordinateRing.getIdeal());
	}

	public Field<T> getField() {
		return field;
	}

	@Override
	public Exactness exactness() {
		return field.exactness();
	}

	@Override
	public AffinePoint<T> getRandomElement() {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean isFinite() {
		return field.isFinite();
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		throw new UnsupportedOperationException();
	}

	@Override
	public Iterator<AffinePoint<T>> iterator() {
		throw new UnsupportedOperationException();
	}

	@Override
	public AffineCover<T> getAffineCover() {
		return cover;
	}

	@Override
	public List<Integer> affineCoverIndex(AffinePoint<T> p) {
		return Collections.singletonList(0);
	}

	public AffineMorphism<T> identityMorphism() {
		List<Polynomial<T>> vars = new ArrayList<>();
		PolynomialRing<T> r = coordinateRing.getPolynomialRing();
		for (int i = 0; i < r.numberOfVariables(); i++) {
			vars.add(r.getVar(i + 1));
		}
		return AffineMorphism.fromPolynomials(this, this, vars);
	}

	@Override
	public Morphism<T, AffinePoint<T>, AffinePoint<T>> embedding(int coverIndex) {
		if (coverIndex != 0) {
			throw new IndexOutOfBoundsException();
		}
		return identityMorphism();
	}

	@Override
	public AffinePoint<T> asAffinePoint(AffinePoint<T> p, int affineCoverIndex) {
		if (affineCoverIndex != 0) {
			throw new IndexOutOfBoundsException();
		}
		return p;
	}

	@Override
	public AffinePoint<T> fromAffinePoint(AffinePoint<T> p, int affineCoverIndex) {
		if (affineCoverIndex != 0) {
			throw new IndexOutOfBoundsException();
		}
		return p;
	}

	@Override
	public Morphism<T, SingletonPoint, AffinePoint<T>> pointAsMorphism(AffinePoint<T> p) {
		SpectrumOfField<T> domain = new SpectrumOfField<>(field);
		return new Morphism<>() {

			@Override
			public AffinePoint<T> evaluate(SingletonPoint t) {
				return p;
			}

			@Override
			public SpectrumOfField<T> getDomain() {
				return domain;
			}

			@Override
			public AffineScheme<T> getRange() {
				return AffineScheme.this;
			}

			@Override
			public RestrictionResult<T> restrict(SingletonPoint preimage) {
				AffineScheme<T> singleton = domain.getAffineCover().getCover().get(0);
				PolynomialRing<T> polynomials = singleton.coordinateRing.getPolynomialRing();
				List<Polynomial<T>> asPolynomials = new ArrayList<>();
				for (int i = 0; i < coordinateRing.getPolynomialRing().numberOfVariables(); i++) {
					asPolynomials.add(polynomials.getEmbedding(p.getCoord(i + 1)));
				}
				return new RestrictionResult<>(0, 0, singleton.identityMorphism(), identityMorphism(),
						AffineMorphism.fromPolynomials(singleton, AffineScheme.this, asPolynomials));
			}

		};
	}
	
	@Override
	public List<AffineScheme<T>> irreducibleComponents() {
		List<AffineScheme<T>> result = new ArrayList<>();
		List<PolynomialIdeal<T>> components = coordinateRing.getIdeal().minimalPrimeIdealsOver();
		for (PolynomialIdeal<T> component : components) {
			result.add(new AffineScheme<>(field, new CoordinateRing<>(coordinateRing.getPolynomialRing(), component)));
		}
		return result;
	}
}
