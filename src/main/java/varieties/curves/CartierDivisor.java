package varieties.curves;

import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import fields.helper.AbstractElement;
import fields.interfaces.Element;
import fields.polynomials.CoordinateRing.CoordinateRingElement;
import varieties.RationalFunction;
import varieties.affine.AffineSubCover;
import varieties.projective.ProjectiveScheme;

//public class CartierDivisor<T extends Element<T>> extends AbstractElement<CartierDivisor<T>> {
//	private ProjectiveScheme<T> scheme;
//	private AffineSubCover<T> cover;
//	private List<RationalFunction<T>> asRationalFunctions;
//
//public LineBundle<T> asLineBundle() {
//	Map<Integer, Map<Integer, CoordinateRingElement<T>>> transitions = new TreeMap<>();
//	for (int i = 0; i < asRationalFunctions.size(); i++) {
//		transitions.put(i, new TreeMap<>());
//		for (int j = 0; j < asRationalFunctions.size(); j++) {
//		RationalFunction<T> transition =	scheme.getFunctionField().divide(asRationalFunctions.get(j), asRationalFunctions.get(i));
//		CoordinateRingElement<T> coord = scheme.getAffineCover().	
//		transitions.get(i).put(j, );
//		}
//	}
//}
//}
