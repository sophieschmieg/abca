package varieties;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import fields.interfaces.Element;
import fields.interfaces.Field;
import util.Pair;
import varieties.GluedScheme.GluedPoint;
import varieties.affine.AffineCover;
import varieties.affine.AffineMorphism;
import varieties.affine.AffineMorphism.AffineFiberedProduct;
import varieties.affine.AffinePoint;
import varieties.affine.AffineScheme;
import varieties.affine.AffineSubCover;
import varieties.affine.AffineSubCover.CommonRefinement;

public class GluedMorphism<T extends Element<T>> extends AbstractMorphism<T, GluedPoint<T>, GluedPoint<T>> {
	private GluedScheme<T> domain;
	private GluedScheme<T> range;
	private AffineSubCover<T> domainSubCover;
	private List<Integer> rangeCoverIndeces;
	private List<List<Integer>> rangeCoverIndexToMorphismIndex;
	private List<AffineMorphism<T>> morphisms;

	public static <T extends Element<T>> GluedMorphism<T> composition(GluedMorphism<T> first, GluedMorphism<T> second) {
		List<List<AffineMorphism<T>>> refinements = new ArrayList<>();
		List<Integer> rangeCoverIndeces = new ArrayList<>();
		List<AffineMorphism<T>> compositions = new ArrayList<>();
		for (int i = 0; i < first.domain.getAffineCover().size(); i++) {
			refinements.add(new ArrayList<>());
		}
		for (int firstFineCoverIndex = 0; firstFineCoverIndex < first.domainSubCover.getFineCover()
				.size(); firstFineCoverIndex++) {
			AffineMorphism<T> morphism = first.morphisms.get(firstFineCoverIndex);
			int intermediateCoarseCoverIndex = first.rangeCoverIndeces.get(firstFineCoverIndex);
			for (int intermediateFineCoverIndex : second.domainSubCover
					.getFineCoverIndex(intermediateCoarseCoverIndex)) {
				AffineMorphism<T> embedding = second.domainSubCover.getRefinement(intermediateFineCoverIndex);
				AffineFiberedProduct<T> product = AffineMorphism.fiberedProduct(embedding, morphism);
				refinements.get(first.domainSubCover.getCoarseCoverIndex(firstFineCoverIndex))
						.add(product.getProjection2());
				compositions.add(AffineMorphism.composition(product.getProjection1(),
						second.morphisms.get(intermediateFineCoverIndex)));
				rangeCoverIndeces.add(second.rangeCoverIndeces.get(intermediateFineCoverIndex));
			}
		}
		return new GluedMorphism<>(first.domain, second.range,
				new AffineSubCover<>(first.domain.getAffineCover(), refinements), rangeCoverIndeces, compositions);
	}

	public GluedMorphism(GluedScheme<T> domain, GluedScheme<T> range, AffineSubCover<T> domainSubCover,
			List<Integer> rangeCoverIndeces, List<AffineMorphism<T>> morphisms) {
		if (domainSubCover.getFineCover().getCover().size() != morphisms.size()
				|| morphisms.size() != rangeCoverIndeces.size()) {
			throw new ArithmeticException("fine cover size mismatch");
		}
		if (domainSubCover.getCoarseCover().getCover().size() != domain.getAffineCover().getCover().size()) {
			throw new ArithmeticException("coarse cover size mismatch");
		}
		for (int i = 0; i < domain.getAffineCover().getCover().size(); i++) {
			if (!domain.getAffineCover().getCover().get(i).equals(domainSubCover.getCoarseCover().getCover().get(i))) {
				throw new ArithmeticException("domain cover mismatch");
			}
		}
		for (int i = 0; i < morphisms.size(); i++) {
			if (!domainSubCover.getFineCover().getCover().get(i).equals(morphisms.get(i).getDomain())) {
				throw new ArithmeticException("domain subcover cover mismatch");
			}
			if (!range.getAffineCover().get(rangeCoverIndeces.get(i)).equals(morphisms.get(i).getRange())) {
				throw new ArithmeticException("range cover mismatch");
			}
		}
		this.domain = domain;
		this.range = range;
		this.domainSubCover = domainSubCover;
		this.rangeCoverIndeces = new ArrayList<>();
		this.rangeCoverIndeces.addAll(rangeCoverIndeces);
		this.morphisms = morphisms;
		this.rangeCoverIndexToMorphismIndex = new ArrayList<>();
		for (int i = 0; i < range.getAffineCover().size(); i++) {
			rangeCoverIndexToMorphismIndex.add(new ArrayList<>());
		}
		for (int i = 0; i < rangeCoverIndeces.size(); i++) {
			rangeCoverIndexToMorphismIndex.get(rangeCoverIndeces.get(i)).add(i);
		}
	}

	public GluedMorphism<T> refine(AffineSubCover<T> refinement) {
		AffineSubCover<T> finerCover = AffineSubCover.composition(refinement, domainSubCover);
		List<AffineMorphism<T>> finerMorphisms = new ArrayList<>();
		List<Integer> finerRangeIndeces = new ArrayList<>();
		for (int finestCoverIndex = 0; finestCoverIndex < refinement.getFineCover().size(); finestCoverIndex++) {
			int fineCoverIndex = refinement.getCoarseCoverIndex(finestCoverIndex);
			finerRangeIndeces.add(rangeCoverIndeces.get(fineCoverIndex));
			finerMorphisms.add(AffineMorphism.composition(refinement.getRefinement(finestCoverIndex),
					morphisms.get(fineCoverIndex)));
		}
		return new GluedMorphism<>(domain, range, finerCover, finerRangeIndeces, finerMorphisms);
	}

	@Override
	public GluedPoint<T> evaluate(GluedPoint<T> t) {
		int domainCoarseCoverIndex = domain.recommendAffineCoverIndex(t);
		AffinePoint<T> preimagePoint = domain.asAffinePoint(t, domainCoarseCoverIndex);
		for (int fineCoverIndex : domainSubCover.getFineCoverIndex(domainCoarseCoverIndex)) {
			List<AffinePoint<T>> refinedPreimagePoints = domainSubCover.getRefinement(fineCoverIndex)
					.preimageList(preimagePoint);
			if (refinedPreimagePoints.isEmpty()) {
				continue;
			}
			if (refinedPreimagePoints.size() > 1) {
				throw new ArithmeticException("Not an injection");
			}
			AffinePoint<T> refinedPreimagePoint = refinedPreimagePoints.get(0);
			AffinePoint<T> imagePoint = morphisms.get(fineCoverIndex).evaluate(refinedPreimagePoint);
			return range.fromAffinePoint(imagePoint, rangeCoverIndeces.get(fineCoverIndex));
		}
		throw new ArithmeticException("Preimage point not found!");
	}

	@Override
	public GluedScheme<T> getDomain() {
		return domain;
	}

	@Override
	public GluedScheme<T> getRange() {
		return range;
	}
	
	public AffineMorphism<T> getMorphism(int i) {
		return morphisms.get(i);
	}
	
	public int subCoverSize() {
		return domainSubCover.getFineCover().size();
	}

	public List<GluedPoint<T>> preimageList(GluedPoint<T> point) {
		Set<GluedPoint<T>> points = new TreeSet<>();
		for (int rangeCoverIndex : range.affineCoverIndex(point)) {
			for (int morphismIndex : rangeCoverIndexToMorphismIndex.get(rangeCoverIndex)) {
				for (AffinePoint<T> preimage : morphisms.get(morphismIndex)
						.preimageList(range.asAffinePoint(point, rangeCoverIndex))) {
					AffinePoint<T> domainPreimage = domainSubCover.getRefinement(morphismIndex).evaluate(preimage);
					points.add(
							domain.fromAffinePoint(domainPreimage, domainSubCover.getCoarseCoverIndex(morphismIndex)));
				}
			}
		}
		List<GluedPoint<T>> result = new ArrayList<>();
		result.addAll(points);
		return result;
	}

	@Override
	public GeneralRationalFunction<T, GluedPoint<T>, GluedPoint<T>> restrict(GluedPoint<T> preimage) {
		int domainCoarseCoverIndex = domain.recommendAffineCoverIndex(preimage);
		AffinePoint<T> preimagePoint = domain.asAffinePoint(preimage, domainCoarseCoverIndex);
		for (int fineCoverIndex : domainSubCover.getFineCoverIndex(domainCoarseCoverIndex)) {
			List<AffinePoint<T>> refinedPreimagePoints = domainSubCover.getRefinement(fineCoverIndex)
					.preimageList(preimagePoint);
			if (refinedPreimagePoints.isEmpty()) {
				continue;
			}
			if (refinedPreimagePoints.size() > 1) {
				throw new ArithmeticException("Not an injection");
			}
			return new GeneralRationalFunction<>(domain, range, morphisms.get(fineCoverIndex), domainCoarseCoverIndex,
					domainSubCover.getRefinement(fineCoverIndex), rangeCoverIndeces.get(fineCoverIndex));
		}
		throw new ArithmeticException("Preimage point not found!");
	}
	
	public GluedMorphism<T> simplifyDomain() {
		return composition(domain.simplify().getInverse(), this);
	}

	public AffineSubCover<T> getDomainSubCover() {
		return domainSubCover;
	}

	public boolean isOpenImmersion() {
		for (AffineMorphism<T> morphism : morphisms) {
			if (!morphism.isOpenImmersion()) {
				return false;
			}
		}
		return true;
	}

	// Proj K[X,Y, Z]/(Y^2Z - X^3 - X^2Z) -> Proj K[X,Y]
	// (X, Z)
	// Spec K[X, Y]/(Y^2-X^3-X^2) --(X)--> Spec K[T]
	// Spec K[X, Z]/(Z-X^3-X^2Z) --(X, Z)--> Spec K[T, S]/(TS-1)  Z(1-X^2)=X^3 X/Z=
	// rangeCoverIndex 0
	// Spec K[X] --(X) -->Spec K[X] x_Spec K[X] id = Spec K[X]
	// Spec K[Y] --(Y)--> Spec K[Y] x_Spec K[X,Y]/(XY-1) Spec K[X] = Spec
	// K[X,Y]/(XY-1)
//	public GluedMorphism<T> restrict(int rangeCoverIndex) {
//		List<AffineScheme<T>> restrictedCover = new ArrayList<>();
//		List<List<AffineMorphism<T>>> restrictedIntersections = new ArrayList<>();
//		List<Integer> restrictedRangeCoverIndeces = new ArrayList<>();
//		List<AffineMorphism<T>> restrictedMorphism = new ArrayList<>();
//		for (int fineDomainCoverIndex = 0; fineDomainCoverIndex < morphisms.size(); fineDomainCoverIndex++) {
//			restrictedRangeCoverIndeces.add(0);
//			AffineFiberedProduct<T> product = AffineMorphism.fiberedProduct(morphisms.get(fineDomainCoverIndex),
//					range.getAffineCover().getIntersectionEmbedding(rangeCoverIndeces.get(fineDomainCoverIndex),
//							rangeCoverIndex));
//			AffineMorphism<T> morphism = AffineMorphism.composition(product.getSimplifiedProjection2(),
//					range.getAffineCover().getIntersectionEmbedding(rangeCoverIndex,
//							rangeCoverIndeces.get(fineDomainCoverIndex)));
//			restrictedCover.add(morphism.getDomain());
//			restrictedMorphism.add(morphism);
//		}
//	}

	public static class GluedFiberedProduct<T extends Element<T>> {
		private GluedScheme<T> product;
		private GluedScheme<T> factor1;
		private GluedScheme<T> factor2;
		private GluedScheme<T> base;
		private GluedMorphism<T> projection1;
		private GluedMorphism<T> projection2;
		private GluedMorphism<T> structural1;
		private GluedMorphism<T> structural2;
		private List<List<AffineFiberedProduct<T>>> affine;
		private List<List<Integer>> productIndeces;
		private List<Pair<Integer, Integer>> factorIndeces;

		public GluedScheme<T> getProduct() {
			return product;
		}

		public GluedScheme<T> getFactor1() {
			return factor1;
		}

		public GluedScheme<T> getFactor2() {
			return factor2;
		}

		public GluedScheme<T> getBase() {
			return base;
		}

		public GluedMorphism<T> getProjection1() {
			return projection1;
		}

		public GluedMorphism<T> getProjection2() {
			return projection2;
		}

		public GluedMorphism<T> getStructural1() {
			return structural1;
		}

		public GluedMorphism<T> getStructural2() {
			return structural2;
		}

		public GluedMorphism<T> universalProperty(GluedMorphism<T> projection1, GluedMorphism<T> projection2) {
			CommonRefinement<T> refinement = new CommonRefinement<>(projection1.domainSubCover,
					projection2.domainSubCover);
			GluedMorphism<T> refinedProjection1 = projection1.refine(refinement.getFirstRefinement());
			GluedMorphism<T> refinedProjection2 = projection2.refine(refinement.getSecondRefinement());
			List<AffineMorphism<T>> morphisms = new ArrayList<>();
			List<Integer> rangeCoverIndeces = new ArrayList<>();
			for (int index = 0; index < refinement.getRefinement().getFineCover().size(); index++) {
				int firstRangeCover = refinedProjection1.rangeCoverIndeces.get(index);
				int secondRangeCover = refinedProjection2.rangeCoverIndeces.get(index);
				rangeCoverIndeces.add(productIndeces.get(firstRangeCover).get(secondRangeCover));
				morphisms.add(affine.get(firstRangeCover).get(secondRangeCover).universalProperty(
						refinedProjection1.morphisms.get(index), refinedProjection2.morphisms.get(index)));
			}
			return new GluedMorphism<>(projection1.getDomain(), getProduct(), refinement.getRefinement(),
					rangeCoverIndeces, morphisms);
		}
	}

	public static <T extends Element<T>> GluedFiberedProduct<T> fiberedProduct(GluedMorphism<T> structural1,
			GluedMorphism<T> structural2) {
		if (structural1.getRange().getAffineCover().size() == 1) {
			return fiberedProductAffineBase(structural1, structural2);
		}
		Field<T> field = structural1.domain.getField();
		GluedFiberedProduct<T> result = new GluedFiberedProduct<>();
		result.structural1 = structural1;
		result.structural2 = structural2;
		result.factor1 = structural1.getDomain();
		result.factor2 = structural2.getDomain();
		result.base = structural1.getRange();
		List<AffineScheme<T>> productCover = new ArrayList<>();
		List<List<AffineMorphism<T>>> productIntersections = new ArrayList<>();
		List<List<List<AffineFiberedProduct<T>>>> affineFiberedProducts = new ArrayList<>();
		List<List<List<Integer>>> productIndeces = new ArrayList<>();
		List<List<Integer>> factorIndeces = new ArrayList<>();
		List<Integer> rangeProjection1CoverIndeces = new ArrayList<>();
		List<AffineMorphism<T>> projection1Morphisms = new ArrayList<>();
		List<Integer> rangeProjection2CoverIndeces = new ArrayList<>();
		List<AffineMorphism<T>> projection2Morphisms = new ArrayList<>();
		AffineCover<T> baseCover = result.base.getAffineCover();
		for (int baseCoverIndex = 0; baseCoverIndex < baseCover.size(); baseCoverIndex++) {
			GluedScheme<T> affineBase = new GluedScheme<>(field, baseCover.get(baseCoverIndex).getAffineCover());
		}
		result.product = new GluedScheme<>(structural1.getDomain().getField(),
				new AffineCover<>(productCover, productIntersections));
		result.projection1 = new GluedMorphism<>(result.product, result.factor1,
				new AffineSubCover<>(result.product.getAffineCover()), rangeProjection1CoverIndeces,
				projection1Morphisms);
		result.projection2 = new GluedMorphism<>(result.product, result.factor2,
				new AffineSubCover<>(result.product.getAffineCover()), rangeProjection2CoverIndeces,
				projection2Morphisms);
//		result.affine = affineFiberedProducts;
//		result.productIndeces = productIndeces;
//		result.factorIndeces = factorIndeces;
		return result;
	}

	private static <T extends Element<T>> GluedFiberedProduct<T> fiberedProductAffineBase(GluedMorphism<T> structural1,
			GluedMorphism<T> structural2) {
		if (structural1.getRange().getAffineCover().size() != 1 || structural2.getRange().getAffineCover().size() != 1
				|| !structural1.getRange().getAffineCover().get(0)
						.equals(structural2.getRange().getAffineCover().get(0))) {
			throw new ArithmeticException("Not an affine base or not the same affine base!");
		}
		GluedFiberedProduct<T> result = new GluedFiberedProduct<>();
		result.structural1 = structural1;
		result.structural2 = structural2;
		result.factor1 = structural1.getDomain();
		result.factor2 = structural2.getDomain();
		result.base = structural1.getRange();
		List<AffineScheme<T>> productCover = new ArrayList<>();
		List<List<AffineMorphism<T>>> productIntersections = new ArrayList<>();
		List<List<AffineFiberedProduct<T>>> affineFiberedProducts = new ArrayList<>();
		List<List<Integer>> productIndeces = new ArrayList<>();
		List<Pair<Integer, Integer>> factorIndeces = new ArrayList<>();
		List<Integer> rangeProjection1CoverIndeces = new ArrayList<>();
		List<AffineMorphism<T>> projection1Morphisms = new ArrayList<>();
		List<Integer> rangeProjection2CoverIndeces = new ArrayList<>();
		List<AffineMorphism<T>> projection2Morphisms = new ArrayList<>();
		int index = 0;
		for (int i = 0; i < structural1.domainSubCover.getFineCover().size(); i++) {
			AffineMorphism<T> structural1Morphism = structural1.morphisms.get(i);
			List<AffineFiberedProduct<T>> affineFiberedProduct = new ArrayList<>();
			List<Integer> productIndex = new ArrayList<>();
			for (int j = 0; j < structural2.domainSubCover.getFineCover().size(); j++) {
				AffineMorphism<T> structural2Morphism = structural2.morphisms.get(j);
				AffineFiberedProduct<T> affineProduct = AffineMorphism.fiberedProduct(structural1Morphism,
						structural2Morphism);
				affineFiberedProduct.add(affineProduct);
				productIndex.add(index);
				factorIndeces.add(new Pair<>(i, j));
				index++;
				productCover.add(affineProduct.getSimplifiedProduct());
				rangeProjection1CoverIndeces.add(structural1.domainSubCover.getCoarseCoverIndex(i));
				projection1Morphisms.add(affineProduct.getSimplifiedProjection1());
				rangeProjection2CoverIndeces.add(structural1.domainSubCover.getCoarseCoverIndex(j));
				projection2Morphisms.add(affineProduct.getSimplifiedProjection2());
			}
			affineFiberedProducts.add(affineFiberedProduct);
			productIndeces.add(productIndex);
		}
		for (int k = 0; k < structural1.domainSubCover.getFineCover().size(); k++) {
			for (int l = 0; l < structural2.domainSubCover.getFineCover().size(); l++) {
				AffineFiberedProduct<T> product = affineFiberedProducts.get(k).get(l);
				List<AffineMorphism<T>> intersections = new ArrayList<>();
				for (int i = 0; i < structural1.domainSubCover.getFineCover().size(); i++) {
					for (int j = 0; j < structural2.domainSubCover.getFineCover().size(); j++) {
						AffineFiberedProduct<T> intersectionProduct = AffineMorphism.fiberedProduct(
								AffineMorphism.composition(
										structural1.domainSubCover.getFineCover().getIntersectionEmbedding(k, i),
										structural1.morphisms.get(k)),
								AffineMorphism.composition(
										structural1.domainSubCover.getFineCover().getIntersectionEmbedding(l, j),
										structural2.morphisms.get(l)));
						AffineMorphism<T> projection1 = AffineMorphism.composition(
								intersectionProduct.getSimplifiedProjection1(),
								structural1.domainSubCover.getFineCover().getIntersectionEmbedding(k, i));
						AffineMorphism<T> projection2 = AffineMorphism.composition(
								intersectionProduct.getSimplifiedProjection2(),
								structural2.domainSubCover.getFineCover().getIntersectionEmbedding(l, j));
						intersections.add(product.simplifiedUniversalProperty(projection1, projection2));
					}
				}
				productIntersections.add(intersections);
			}
		}
		result.product = new GluedScheme<>(structural1.getDomain().getField(),
				new AffineCover<>(productCover, productIntersections));
		result.projection1 = new GluedMorphism<>(result.product, result.factor1,
				new AffineSubCover<>(result.product.getAffineCover()), rangeProjection1CoverIndeces,
				projection1Morphisms);
		result.projection2 = new GluedMorphism<>(result.product, result.factor2,
				new AffineSubCover<>(result.product.getAffineCover()), rangeProjection2CoverIndeces,
				projection2Morphisms);
		result.affine = affineFiberedProducts;
		result.productIndeces = productIndeces;
		result.factorIndeces = factorIndeces;
		return result;
	}

}
