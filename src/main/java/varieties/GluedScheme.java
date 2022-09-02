package varieties;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Deque;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Optional;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;

import fields.exceptions.InfinityException;
import fields.helper.AbstractElement;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import varieties.GluedScheme.GluedPoint;
import varieties.SpectrumOfField.SingletonPoint;
import varieties.affine.AffineCover;
import varieties.affine.AffineMorphism;
import varieties.affine.AffineMorphism.AffineFiberedProduct;
import varieties.affine.AffinePoint;
import varieties.affine.AffineScheme;
import varieties.affine.AffineScheme.AffineNormalization;
import varieties.affine.AffineScheme.AffineSimplification;
import varieties.affine.AffineSubCover;
import varieties.affine.AffineSubCover.CommonRefinement;

public class GluedScheme<T extends Element<T>> extends AbstractScheme<T, GluedPoint<T>>
		implements Scheme<T, GluedPoint<T>> {
	public static class GluedPoint<T extends Element<T>> extends AbstractElement<GluedPoint<T>> {
		private AffinePoint<T> point;
		private int affineCoverIndex;

		private GluedPoint(AffinePoint<T> point, int affineCoverIndex) {
			this.point = point;
			this.affineCoverIndex = affineCoverIndex;
		}

		@Override
		public int compareTo(GluedPoint<T> o) {
			if (affineCoverIndex != o.affineCoverIndex) {
				return o.affineCoverIndex - affineCoverIndex;
			}
			return point.compareTo(o.point);
		}

		@Override
		public String toString() {
			return point + "@" + affineCoverIndex;
		}
	}

	public static <T extends Element<T>> GluedScheme<T> glue(AffineMorphism<T> t1, AffineMorphism<T> t2) {
		if (!t1.isOpenImmersion() || !t2.isOpenImmersion()) {
			throw new ArithmeticException("Not open immersions!");
		}
		if (!t1.getDomain().equals(t2.getDomain())) {
			throw new ArithmeticException("No shared domain!");
		}
		List<AffineScheme<T>> cover = new ArrayList<>();
		List<List<AffineMorphism<T>>> intersections = new ArrayList<>();
		cover.add(t1.getRange());
		cover.add(t2.getRange());
		intersections.add(new ArrayList<>());
		intersections.get(0).add(t1.getRange().identityMorphism());
		intersections.get(0).add(t1);
		intersections.add(new ArrayList<>());
		intersections.get(1).add(t2);
		intersections.get(1).add(t2.getRange().identityMorphism());
		return new GluedScheme<>(t1.getDomain().getField(), new AffineCover<>(cover, intersections));
	}

	public static <T extends Element<T>> GluedScheme<T> glue(GluedMorphism<T> t1, GluedMorphism<T> t2) {
		if (!t1.isOpenImmersion() || !t2.isOpenImmersion()) {
			throw new ArithmeticException("Not open immersions!");
		}
		if (!t1.getDomain().equals(t2.getDomain())) {
			throw new ArithmeticException("No shared domain!");
		}
		CommonRefinement<T> refinement =new AffineSubCover.CommonRefinement<>(t1.getDomainSubCover(), t2.getDomainSubCover());
		t1 = t1.refine(refinement.getFirstRefinement());
		t2 = t2.refine(refinement.getSecondRefinement());
		List<AffineScheme<T>> cover = new ArrayList<>();
		List<List<AffineMorphism<T>>> intersections = new ArrayList<>();
		for (int i = 0; i < t1.subCoverSize(); i++) {
			cover.add(t1.getMorphism(i).getRange());
			List<AffineMorphism<T>> intersection = new ArrayList<>();
			for (int j = 0; j < t1.subCoverSize(); j++) {
				intersection.add(t1.getRange().getAffineCover().getIntersectionEmbedding(i, j));
			}
			for (int j = 0; j < t2.subCoverSize(); j++) {
			}
			intersections.add(intersection);
		}
		for (int i = 0; i < t2.subCoverSize(); i++) {
			cover.add(t2.getMorphism(i).getRange());
		}
		return null;
	}

	private Field<T> field;
	private AffineCover<T> affineCover;
	private GluedSimplification<T> simplification;

	public GluedScheme(Field<T> field, AffineCover<T> affineCover) {
		this.field = field;
		this.affineCover = affineCover;
	}

	@Override
	public Exactness exactness() {
		return field.exactness();
	}

	private Optional<AffinePoint<T>> getPoint(AffinePoint<T> point, int affineCoverIndex, int targetCoverIndex) {
		if (targetCoverIndex == affineCoverIndex) {
			return Optional.of(point);
		}
		AffineMorphism<T> intersectionToAffineCoverIndex = affineCover.getIntersectionEmbedding(affineCoverIndex,
				targetCoverIndex);
		AffineMorphism<T> intersectionToCover = affineCover.getIntersectionEmbedding(targetCoverIndex,
				affineCoverIndex);
		List<AffinePoint<T>> points = new ArrayList<>();
		for (AffinePoint<T> preImagePoint : intersectionToAffineCoverIndex.preimageList(point)) {
			points.add(intersectionToCover.evaluate(preImagePoint));
		}
		if (points.isEmpty()) {
			return Optional.empty();
		}
		if (points.size() != 1) {
			throw new ArithmeticException("more than one point in preimage!");
		}
		return Optional.of(points.get(0));
	}

	@Override
	public GluedPoint<T> getRandomElement() {
		int coverIndex = new Random().nextInt(affineCover.getCover().size());
		AffinePoint<T> point = affineCover.getCover().get(coverIndex).getRandomElement();
		return fromAffinePoint(point, coverIndex);
	}

	@Override
	public boolean isFinite() {
		return field.isFinite() || dimension() == 0;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		throw new UnsupportedOperationException("Point counting not possible!");
	}

	@Override
	public Iterator<GluedPoint<T>> iterator() {
		return new Iterator<>() {
			int affineCoverIndex = affineCover.getCover().size() - 1;
			Iterator<AffinePoint<T>> it = affineCover.getCover().get(affineCoverIndex).iterator();
			Deque<GluedPoint<T>> queue = new LinkedList<>();

			private boolean fillQueue() {
				while (queue.isEmpty()) {
					while (!it.hasNext()) {
						if (affineCoverIndex == 0) {
							return false;
						}
						affineCoverIndex--;
						it = affineCover.getCover().get(affineCoverIndex).iterator();
					}
					AffinePoint<T> point = it.next();
					GluedPoint<T> gluedPoint = fromAffinePoint(point, affineCoverIndex);
					if (gluedPoint.affineCoverIndex == affineCoverIndex) {
						queue.add(gluedPoint);
					}
				}
				return true;
			}

			@Override
			public boolean hasNext() {
				return fillQueue();
			}

			@Override
			public GluedPoint<T> next() {
				fillQueue();
				return queue.poll();
			}
		};
	}

	@Override
	public Field<T> getField() {
		return field;
	}

	@Override
	public boolean hasRationalPoint(GluedPoint<T> p) {
		return affineCover.getCover().get(p.affineCoverIndex).hasRationalPoint(p.point);
	}

	@Override
	public AffineCover<T> getAffineCover() {
		return affineCover;
	}

	@Override
	public int recommendAffineCoverIndex(GluedPoint<T> p) {
		return p.affineCoverIndex;
	}

	@Override
	public List<Integer> affineCoverIndex(GluedPoint<T> p) {
		List<Integer> result = new ArrayList<>();
		for (int i = affineCover.getCover().size() - 1; i >= 0; i--) {
			if (getPoint(p.point, p.affineCoverIndex, i).isPresent()) {
				result.add(i);
			}
		}
		return result;
	}

	@Override
	public AffinePoint<T> asAffinePoint(GluedPoint<T> p, int affineCoverIndex) {
		return getPoint(p.point, p.affineCoverIndex, affineCoverIndex).get();
	}

	@Override
	public GluedPoint<T> fromAffinePoint(AffinePoint<T> p, int affineCoverIndex) {
		for (int i = affineCover.getCover().size() - 1; i >= 0; i--) {
			Optional<AffinePoint<T>> converted = getPoint(p, affineCoverIndex, i);
			if (converted.isPresent()) {
				return new GluedPoint<>(converted.get(), i);
			}
		}
		throw new ArithmeticException("Point not found!");
	}

	@Override
	public GluedMorphism<T> identityMorphism() {
		List<AffineMorphism<T>> identities = new ArrayList<>();
		List<Integer> rangeCoverIndeces = new ArrayList<>();
		for (int i = 0; i < affineCover.size(); i++) {
			identities.add(affineCover.get(i).identityMorphism());
			rangeCoverIndeces.add(i);
		}
		return new GluedMorphism<>(this, this, new AffineSubCover<>(affineCover), rangeCoverIndeces, identities);
	}

	@Override
	public Morphism<T, AffinePoint<T>, GluedPoint<T>> embedding(int affineCoverIndex) {
		return new AbstractMorphism<>() {

			@Override
			public GluedPoint<T> evaluate(AffinePoint<T> t) {
				return fromAffinePoint(t, affineCoverIndex);
			}

			@Override
			public Scheme<T, AffinePoint<T>> getDomain() {
				return affineCover.getCover().get(affineCoverIndex);
			}

			@Override
			public GluedScheme<T> getRange() {
				return GluedScheme.this;
			}

			@Override
			public GeneralRationalFunction<T, AffinePoint<T>, GluedPoint<T>> restrict(AffinePoint<T> point) {
				return new GeneralRationalFunction<>(affineCover.getCover().get(affineCoverIndex), GluedScheme.this,
						affineCover.get(affineCoverIndex).identityMorphism(), 0,
						affineCover.get(affineCoverIndex).identityMorphism(), affineCoverIndex);
			}
		};
	}

	@Override
	public Morphism<T, SingletonPoint, GluedPoint<T>> pointAsMorphism(GluedPoint<T> p) {
		SpectrumOfField<T> domain = new SpectrumOfField<>(field);
		return new AbstractMorphism<>() {

			@Override
			public GluedPoint<T> evaluate(SingletonPoint t) {
				return p;
			}

			@Override
			public SpectrumOfField<T> getDomain() {
				return domain;
			}

			@Override
			public GluedScheme<T> getRange() {
				return GluedScheme.this;
			}

			@Override
			public GeneralRationalFunction<T, SingletonPoint, GluedPoint<T>> restrict(SingletonPoint point) {
				int rangeCoverIndex = recommendAffineCoverIndex(p);
				AffineScheme<T> singleton = domain.getAffineCover().getCover().get(0);
				PolynomialRing<T> polynomials = singleton.getCoordinateRing().getPolynomialRing();
				List<Polynomial<T>> asPolynomials = new ArrayList<>();
				AffineScheme<T> affineRange = getAffineCover().getCover().get(rangeCoverIndex);
				for (T coord : asAffinePoint(p, rangeCoverIndex).getCoords()) {
					asPolynomials.add(polynomials.getEmbedding(coord));
				}
				return new GeneralRationalFunction<>(domain, GluedScheme.this,
						AffineMorphism.fromPolynomials(singleton, affineRange, asPolynomials), 0,
						singleton.identityMorphism(), rangeCoverIndex);
			}

		};
	}

	public static class GluedSimplification<T extends Element<T>> {
		private GluedMorphism<T> simplification;
		private GluedMorphism<T> inverse;

		private GluedSimplification(GluedMorphism<T> simplification, GluedMorphism<T> inverse) {
			this.simplification = simplification;
			this.inverse = inverse;
		}

		public GluedMorphism<T> getSimplification() {
			return simplification;
		}

		public GluedMorphism<T> getInverse() {
			return inverse;
		}
	}

	private GluedSimplification<T> simplifyCover() {
		List<AffineScheme<T>> simplerCover = new ArrayList<>();
		List<List<AffineMorphism<T>>> simplerIntersections = new ArrayList<>();
		List<Integer> simplerRangeIndeces = new ArrayList<>();
		List<AffineMorphism<T>> simplerMorphism = new ArrayList<>();
		List<AffineMorphism<T>> simplerInverseMorphism = new ArrayList<>();
		AffineCover<T> simplerIntersectionsCover = affineCover.simplifyIntersections();
		for (int i = 0; i < affineCover.size(); i++) {
			AffineSimplification<T> simplification = affineCover.get(i).simplify();
			simplerCover.add(simplification.getSimplification().getRange());
			simplerMorphism.add(simplification.getSimplification());
			simplerInverseMorphism.add(simplification.getInverse());
			simplerIntersections.add(new ArrayList<>());
			simplerRangeIndeces.add(i);
			for (int j = 0; j < affineCover.size(); j++) {
				simplerIntersections.get(i).add(AffineMorphism.composition(
						simplerIntersectionsCover.getIntersectionEmbedding(i, j), simplification.getSimplification()));
			}
		}
		GluedScheme<T> simplified = new GluedScheme<>(field, new AffineCover<>(simplerCover, simplerIntersections));
		GluedMorphism<T> simplfication = new GluedMorphism<>(this, simplified, new AffineSubCover<>(affineCover),
				simplerRangeIndeces, simplerMorphism);
		GluedMorphism<T> inverse = new GluedMorphism<>(simplified, this,
				new AffineSubCover<>(simplified.getAffineCover()), simplerRangeIndeces, simplerInverseMorphism);
		return new GluedSimplification<>(simplfication, inverse);
	}

	private GluedSimplification<T> simplifySubsets() {
		List<AffineScheme<T>> simplerCover = new ArrayList<>();
		List<List<AffineMorphism<T>>> simplerIntersections = new ArrayList<>();
		List<Integer> simplerRangeIndeces = new ArrayList<>();
		List<AffineMorphism<T>> simplerMorphism = new ArrayList<>();
		List<Integer> simplerInverseRangeIndeces = new ArrayList<>();
		List<AffineMorphism<T>> simplerInverseMorphism = new ArrayList<>();
		List<List<Integer>> subsets = new ArrayList<>();
		for (int i = 0; i < affineCover.size(); i++) {
			simplerRangeIndeces.add(i);
			simplerMorphism.add(affineCover.get(i).identityMorphism());
			List<Integer> subset = new ArrayList<>();
			subset.add(i);
			subsets.add(subset);
		}
		Set<Integer> subsetIndex = new TreeSet<>();
		for (int i = 0; i < affineCover.size(); i++) {
			// U --> U ^ V --> V
			boolean redundant = false;
			for (int j = 0; j < affineCover.size(); j++) {
				if (i == j || subsetIndex.contains(j)) {
					continue;
				}
				Optional<AffineMorphism<T>> inverse = affineCover.getIntersectionEmbedding(i, j).inverse();
				if (inverse.isPresent()) {
					for (int k : subsets.get(i)) {
						simplerMorphism.set(k, AffineMorphism.composition(simplerMorphism.get(k),
								AffineMorphism.composition(inverse.get(), affineCover.getIntersectionEmbedding(j, i))));
						simplerRangeIndeces.set(k, j);
						subsets.get(j).add(k);
					}
					subsetIndex.add(i);
					redundant = true;
					break;
				}
			}
			if (!redundant) {
				int simplifiedIndex = simplerCover.size();
				simplerCover.add(affineCover.get(i));
				simplerInverseRangeIndeces.add(i);
				simplerInverseMorphism.add(affineCover.get(i).identityMorphism());
				for (int j : subsets.get(i)) {
					simplerRangeIndeces.set(j, simplifiedIndex);
				}
			}
		}
		for (int i = 0; i < affineCover.size(); i++) {
			if (subsetIndex.contains(i)) {
				continue;
			}
			List<AffineMorphism<T>> intersection = new ArrayList<>();
			for (int j = 0; j < affineCover.size(); j++) {
				if (subsetIndex.contains(j)) {
					continue;
				}
				intersection.add(affineCover.getIntersectionEmbedding(i, j));
			}
			simplerIntersections.add(intersection);
		}
		GluedScheme<T> simplified = new GluedScheme<>(field, new AffineCover<>(simplerCover, simplerIntersections));
		GluedMorphism<T> simplification = new GluedMorphism<>(this, simplified, new AffineSubCover<>(affineCover),
				simplerRangeIndeces, simplerMorphism);
		GluedMorphism<T> inverse = new GluedMorphism<>(simplified, this,
				new AffineSubCover<>(simplified.getAffineCover()), simplerInverseRangeIndeces, simplerInverseMorphism);
		return new GluedSimplification<>(simplification, inverse);
	}

	public GluedSimplification<T> simplify() {
		if (simplification == null) {
			GluedSimplification<T> subsetSimplification = simplifySubsets();
			GluedSimplification<T> coverSimplificiation = subsetSimplification.getSimplification().getRange()
					.simplifyCover();
			simplification = new GluedSimplification<>(
					GluedMorphism.composition(subsetSimplification.getSimplification(),
							coverSimplificiation.getSimplification()),
					GluedMorphism.composition(coverSimplificiation.getInverse(), subsetSimplification.getInverse()));
		}
		return simplification;
	}

	public GluedMorphism<T> normalization() {
		GluedSimplification<T> simpler = simplify();
		return GluedMorphism
				.composition(simpler.getSimplification().getRange().computeNormalization(), simpler.getInverse())
				.simplifyDomain();
	}

	public GluedMorphism<T> computeNormalization() {
		List<AffineScheme<T>> affineNormalizationDomains = new ArrayList<>();
		List<AffineMorphism<T>> restrictions = new ArrayList<>();
		List<AffineMorphism<T>> affineNormalizations = new ArrayList<>();
		List<AffineNormalization<T>> normalizations = new ArrayList<>();
		List<Integer> rangeIndeces = new ArrayList<>();
		List<List<AffineMorphism<T>>> affineNormalizationAtlas = new ArrayList<>();
		for (int index = 0; index < affineCover.size(); index++) {
			AffineScheme<T> affine = affineCover.get(index);
			AffineNormalization<T> normalization = affine.normalization();
			normalizations.add(normalization);
			affineNormalizations.add(normalization.getNormalizationMorphism());
			rangeIndeces.add(index);
			affineNormalizationDomains.add(normalization.getNormalization());
			AffineMorphism<T> inverseDomainEmbedding = normalization.getInverse().getDomainEmbedding();
			restrictions.add(AffineMorphism.getClosedImmersion(inverseDomainEmbedding.getRange(),
					inverseDomainEmbedding.getRange().getCoordinateRing()
							.getIdeal(Collections.singletonList(inverseDomainEmbedding.getInversionObstruction()))));
			List<AffineMorphism<T>> intersections = new ArrayList<>();
			for (int i = 0; i < affineCover.size(); i++) {
				intersections.add(null);
			}
			affineNormalizationAtlas.add(intersections);
			affineNormalizationAtlas.get(index).set(index, affineNormalizationDomains.get(index).identityMorphism());
		}
		for (int i = 1; i < affineCover.size(); i++) {
			for (int j = 0; j < i; j++) {
				AffineMorphism<T> intersectionToI = affineCover.getIntersectionEmbedding(i, j);
				AffineMorphism<T> intersectionToJ = affineCover.getIntersectionEmbedding(j, i);
				AffineFiberedProduct<T> intersectionToNormalizationI = AffineMorphism.fiberedProduct(intersectionToI,
						affineNormalizations.get(i));
				AffineFiberedProduct<T> intersectionToNormalizationJ = AffineMorphism.fiberedProduct(intersectionToJ,
						affineNormalizations.get(j));
				AffineFiberedProduct<T> intersection = AffineMorphism.fiberedProduct(
						intersectionToNormalizationI.getSimplifiedProjection1(),
						intersectionToNormalizationJ.getSimplifiedProjection1());
				affineNormalizationAtlas.get(i).set(j,
						AffineMorphism.composition(intersection.getSimplifiedProjection1(),
								intersectionToNormalizationI.getSimplifiedProjection2()));
				affineNormalizationAtlas.get(j).set(i,
						AffineMorphism.composition(intersection.getSimplifiedProjection2(),
								intersectionToNormalizationJ.getSimplifiedProjection2()));
			}
		}
		AffineCover<T> normalizationCover = new AffineCover<>(affineNormalizationDomains, affineNormalizationAtlas);
		GluedScheme<T> normalization = new GluedScheme<>(field, normalizationCover);
		GluedMorphism<T> normalizationMorphism = new GluedMorphism<>(normalization, this,
				new AffineSubCover<>(normalization.getAffineCover()), rangeIndeces, affineNormalizations);
		return normalizationMorphism;
	}

	@Override
	public List<GluedMorphism<T>> irreducibleComponents() {
		List<GluedMorphism<T>> result = new ArrayList<>();
		List<List<AffineMorphism<T>>> irreducibleComponents = new ArrayList<>();
		for (int i = 0; i < affineCover.size(); i++) {
			irreducibleComponents.add(affineCover.get(i).irreducibleComponents());
		}
		for (int i = 0; i < affineCover.size(); i++) {
			List<AffineMorphism<T>> components = irreducibleComponents.get(i);
			componentLoop: for (AffineMorphism<T> component : components) {
				for (int j = 0; j < i; j++) {
					AffineMorphism<T> inIntersection = affineCover.getIntersectionEmbedding(i, j).preimage(component);
					if (inIntersection.getDomain().dimension() >= 0) {
						continue componentLoop;
					}
				}
				List<AffineMorphism<T>> componentMorphisms = new ArrayList<>();
				List<AffineScheme<T>> componentSlices = new ArrayList<>();
				List<List<AffineMorphism<T>>> componentIntersections = new ArrayList<>();
				for (int j = 0; j < affineCover.size(); j++) {
					componentIntersections.add(new ArrayList<>());
					for (int k = 0; k < affineCover.size(); k++) {
						componentIntersections.get(j).add(null);
					}
				}
				for (int j = 0; j < affineCover.size(); j++) {
					AffineMorphism<T> intersectionToI = affineCover.getIntersectionEmbedding(i, j);
					AffineMorphism<T> componentToIntersection = intersectionToI.preimage(component);
					AffineMorphism<T> intersectionToJ = affineCover.getIntersectionEmbedding(j, i);
					AffineMorphism<T> iComponentToJ = AffineMorphism.composition(componentToIntersection,
							intersectionToJ);
					AffineMorphism<T> jComponent = iComponentToJ.closure();
					componentSlices.add(jComponent.getDomain());
					componentMorphisms.add(jComponent);
					for (int k = 0; k < affineCover.size(); k++) {
						AffineMorphism<T> intersectionJWithK = affineCover.getIntersectionEmbedding(j, k);
						componentIntersections.get(j).set(k, jComponent.preimage(intersectionJWithK));
					}
				}
				GluedScheme<T> componentScheme = new GluedScheme<>(field,
						new AffineCover<>(componentSlices, componentIntersections));
//				result.add(new GluedMorphism<>(componentScheme, this,
//						new AffineSubCover<>(componentScheme.getAffineCover()), new AffineSubCover<>(affineCover),
//						componentMorphisms));
			}
		}
		return result;
	}

	@Override
	public GluedMorphism<T> reduced() {
		// TODO Auto-generated method stub
		return null;
	}
}
