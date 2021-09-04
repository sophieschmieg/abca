package util;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

class CacheTest {

	@Test
	void test() {
		Cache<Integer, String> cache = new Cache<>(5);
		assertEquals(0, cache.size());
		cache.insert(1, "A");
		assertEquals(1, cache.size());
		assertEquals("A", cache.lookup(1).get());
		assertEquals(1, cache.size());
		assertEquals("A", cache.lookup(1).get());
		assertTrue(cache.lookup(2).isEmpty());
		assertEquals(1, cache.size());
		cache.remove(1);
		assertEquals(0, cache.size());
		assertTrue(cache.lookup(1).isEmpty());
		cache.insert(3, "C");
		cache.insert(4, "D");
		assertEquals(2, cache.size());
		assertEquals("C", cache.lookup(3).get());
		assertEquals("D", cache.lookup(4).get());
		cache.remove(3);
		cache.remove(4);
		assertEquals(0, cache.size());
		assertTrue(cache.lookup(3).isEmpty());
		assertTrue(cache.lookup(4).isEmpty());
		cache.insert(1, "a");
		cache.insert(2, "b");
		cache.insert(3, "c");
		cache.insert(4, "d");
		cache.insert(5, "e");
		cache.insert(6, "f");
		assertEquals(5, cache.size());
		assertTrue(cache.lookup(1).isEmpty());
		assertEquals("b", cache.lookup(2).get());
		cache.insert(7, "g");
		assertEquals("b", cache.lookup(2).get());
		assertTrue(cache.lookup(3).isEmpty());
	}

}
