package varieties.projective;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.CoordinateRing;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import varieties.AbstractScheme;
import varieties.Morphism;
import varieties.SpectrumOfField;
import varieties.SpectrumOfField.SingletonPoint;
import varieties.affine.AffineCover;
import varieties.affine.AffineMorphism;
import varieties.affine.AffinePoint;
import varieties.affine.AffineScheme;

public class ProjectiveScheme<T extends Element<T>> extends AbstractScheme<T, ProjectivePoint<T>> {
	private PolynomialRing<T> polynomialRing;
	private List<Polynomial<T>> generators;
	private Field<T> field;
	private AffineCover<T> cover;

	public ProjectiveScheme(Field<T> field, PolynomialRing<T> polynomialRing, List<Polynomial<T>> generators) {
		this.field = field;
		this.polynomialRing = polynomialRing;
		this.generators = polynomialRing.getIdeal(generators).generators();
		for (Polynomial<T> generator : generators) {
			if (!polynomialRing.isHomogeneous(generator)) {
				throw new ArithmeticException("not homogenous");
			}
		}
		List<AffineScheme<T>> coverList = new ArrayList<>();
		List<List<AffineMorphism<T>>> intersectionList = new ArrayList<>();
		PolynomialRing<T> affineRing = AbstractPolynomialRing.getPolynomialRing(field,
				polynomialRing.numberOfVariables() - 1, polynomialRing.getComparator());
		Polynomial<T> z = polynomialRing.getVar(polynomialRing.numberOfVariables());
		for (int i = 0; i < polynomialRing.numberOfVariables(); i++) {
			List<Polynomial<T>> affineList = new ArrayList<>();
			int[] map = new int[polynomialRing.numberOfVariables()];
			for (int j = 0; j < polynomialRing.numberOfVariables(); j++) {
				if (j < i) {
					map[j] = j;
				} else if (j > i) {
					map[j] = j - 1;
				} else {
					map[j] = -1;
				}
			}
			for (Polynomial<T> generator : generators) {
				affineList.add(affineRing.getEmbedding(polynomialRing.dehomogenize(generator, i + 1), map));
			}
			CoordinateRing<T> cr = new CoordinateRing<>(affineRing, affineRing.getIdeal(affineList));
			coverList.add(new AffineScheme<>(field, cr));
			intersectionList.add(new ArrayList<>());
			for (int j = 0; j < i; j++) {
				Polynomial<T> inverse = polynomialRing
						.subtract(polynomialRing.multiply(polynomialRing.getVar(j + 1), z), polynomialRing.one());
				List<Polynomial<T>> intersectionGenerators = new ArrayList<>();
				intersectionGenerators.add(inverse);
				for (Polynomial<T> generator : affineList) {
					intersectionGenerators.add(polynomialRing.getEmbedding(generator));
				}
				CoordinateRing<T> intersectionCoordinateRing = new CoordinateRing<>(polynomialRing,
						polynomialRing.getIdeal(intersectionGenerators));
				AffineScheme<T> intersection = new AffineScheme<>(field, intersectionCoordinateRing);
				List<Polynomial<T>> firstEmbeddingList = new ArrayList<>();
				List<Polynomial<T>> secondEmbeddingList = new ArrayList<>();
				for (int k = 0; k < affineRing.numberOfVariables(); k++) {
					firstEmbeddingList.add(polynomialRing.getVar(k + 1));
					Polynomial<T> var;
					if (k < j || k >= i) {
						var = polynomialRing.getVar(k + 1);
					} else if (k == i - 1) {
						var = polynomialRing.one();
					} else {
						var = polynomialRing.getVar(k + 2);
					}
					secondEmbeddingList.add(polynomialRing.multiply(var, z));
				}
				AffineMorphism<T> firstEmbedding = AffineMorphism.fromPolynomials(intersection, coverList.get(i),
						firstEmbeddingList);
				AffineMorphism<T> secondEmbedding = AffineMorphism.fromPolynomials(intersection, coverList.get(j),
						secondEmbeddingList);
				intersectionList.get(i).add(firstEmbedding);
				intersectionList.get(j).add(secondEmbedding);
			}
			intersectionList.get(i).add(coverList.get(i).identityMorphism());
		}
		this.cover = new AffineCover<>(coverList, intersectionList);
	}

	public static <T extends Element<T>> ProjectiveScheme<T> fromAffineScheme(AffineScheme<T> affine) {
		PolynomialRing<T> affinePolynomialRing = affine.getCoordinateRing().getPolynomialRing();
		PolynomialRing<T> homogenousPolynomialRing = AbstractPolynomialRing.getPolynomialRing(affine.getField(),
				affinePolynomialRing.numberOfVariables() + 1, affinePolynomialRing.getComparator());
		List<Polynomial<T>> polynomials = new ArrayList<>();
		for (Polynomial<T> polynominal : affine.getCoordinateRing().getIdeal().generators()) {
			polynomials.add(affinePolynomialRing.homogenize(polynominal));
		}
		return new ProjectiveScheme<>(affine.getField(), homogenousPolynomialRing, polynomials);
	}

	public PolynomialRing<T> homogenousPolynomialRing() {
		return polynomialRing;
	}

	@Override
	public Exactness exactness() {
		return field.exactness();
	}

	@Override
	public ProjectivePoint<T> getRandomElement() {
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
	public Iterator<ProjectivePoint<T>> iterator() {
		throw new UnsupportedOperationException();
	}

	@Override
	public Field<T> getField() {
		return field;
	}

	@Override
	public boolean hasRationalPoint(ProjectivePoint<T> p) {
		for (Polynomial<T> generator : generators) {
			if (!polynomialRing.evaluate(generator, p.getCoords()).equals(field.zero())) {
				return false;
			}
		}
		return true;
	}

	@Override
	public int dimension() {
		return cover.getCover().get(0).dimension();
	}

	public int projectiveEmbeddingDimension() {
		return polynomialRing.numberOfVariables() - 1;
	}

	@Override
	public AffineCover<T> getAffineCover() {
		return cover;
	}

	@Override
	public List<Integer> affineCoverIndex(ProjectivePoint<T> p) {
		List<Integer> result = new ArrayList<>();
		for (int i = polynomialRing.numberOfVariables() - 1; i >= 0; i--) {
			T coord = p.getCoord(i + 1);
			if (!coord.equals(field.zero())) {
				result.add(i);
			}
		}
		return result;
	}

	@Override
	public Morphism<T, AffinePoint<T>, ProjectivePoint<T>> embedding(int coverIndex) {
		return new Morphism<>() {

			@Override
			public ProjectivePoint<T> evaluate(AffinePoint<T> t) {
				List<T> values = new ArrayList<>();
				for (int i = 0; i < polynomialRing.numberOfVariables(); i++) {
					if (i == coverIndex) {
						values.add(field.one());
					} else {
						values.add(t.getCoord(i < coverIndex ? i + 1 : i));
					}
				}
				return new ProjectivePoint<>(field, values);
			}

			@Override
			public AffineScheme<T> getDomain() {
				return cover.getCover().get(coverIndex);
			}

			@Override
			public ProjectiveScheme<T> getRange() {
				return ProjectiveScheme.this;
			}

			@Override
			public RestrictionResult<T> restrict(AffinePoint<T> p) {
				return new RestrictionResult<>(coverIndex, coverIndex,
						cover.getCover().get(coverIndex).identityMorphism(),
						cover.getCover().get(coverIndex).identityMorphism(),
						cover.getCover().get(coverIndex).identityMorphism());
			}
		};
	}

	@Override
	public AffinePoint<T> asAffinePoint(ProjectivePoint<T> p, int affineCoverIndex) {
		return p.getDehomogenous(affineCoverIndex + 1);
	}

	@Override
	public ProjectivePoint<T> fromAffinePoint(AffinePoint<T> p, int affineCoverIndex) {
		List<T> coords = new ArrayList<>();
		for (int i = 0; i < polynomialRing.numberOfVariables(); i++) {
			if (i == affineCoverIndex) {
				coords.add(field.one());
				continue;
			}
			coords.add(p.getCoord(i > affineCoverIndex ? i : i + 1));
		}
		return new ProjectivePoint<>(field, coords);
	}

	@Override
	public Morphism<T, SingletonPoint, ProjectivePoint<T>> pointAsMorphism(ProjectivePoint<T> p) {
		SpectrumOfField<T> domain = new SpectrumOfField<>(field);
		return new Morphism<>() {

			@Override
			public ProjectivePoint<T> evaluate(SingletonPoint t) {
				return p;
			}

			@Override
			public SpectrumOfField<T> getDomain() {
				return domain;
			}

			@Override
			public ProjectiveScheme<T> getRange() {
				return ProjectiveScheme.this;
			}

			@Override
			public RestrictionResult<T> restrict(SingletonPoint preimage) {
				AffineScheme<T> singleton = domain.getAffineCover().getCover().get(0);
				PolynomialRing<T> polynomials = singleton.getCoordinateRing().getPolynomialRing();
				List<Polynomial<T>> asPolynomials = new ArrayList<>();
				int rangeCoverIndex = affineCoverIndex(p).get(0);
				AffineScheme<T> affineRange = getAffineCover().getCover().get(rangeCoverIndex);
				for (int i = 0; i < polynomialRing.numberOfVariables(); i++) {
					if (i == rangeCoverIndex) {
						continue;
					}
					asPolynomials.add(polynomials.getEmbedding(p.getDehomogenisedCoord(i + 1, rangeCoverIndex + 1)));
				}
				return new RestrictionResult<>(0, rangeCoverIndex, singleton.identityMorphism(),
						affineRange.identityMorphism(),
						AffineMorphism.fromPolynomials(singleton, affineRange, asPolynomials));
			}

		};
	}

	@Override
	public List<ProjectiveScheme<T>> irreducibleComponents() {
		List<ProjectiveScheme<T>> result = new ArrayList<>();
		for (AffineScheme<T> affine : getAffineCover().getCover().get(0).irreducibleComponents()) {
			result.add(fromAffineScheme(affine));
		}
		return result;
	}
}
