﻿-------------------------------------------------------------------------------

DROP INDEX IF EXISTS tbl_tgeompoint3D_spgist_idx;
DROP INDEX IF EXISTS tbl_tgeogpoint3D_spgist_idx;

-------------------------------------------------------------------------------

ALTER TABLE test_geoboundboxops3d ADD spgistidx bigint ;

-------------------------------------------------------------------------------

CREATE INDEX tbl_tgeompoint3D_spgist_idx ON tbl_tgeompoint3D USING SPGIST(temp);
CREATE INDEX tbl_tgeogpoint3D_spgist_idx ON tbl_tgeogpoint3D USING SPGIST(temp);

-------------------------------------------------------------------------------
-- <type> op tgeompoint3D

UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_geometry3D, tbl_tgeompoint3D WHERE g && temp )
WHERE op = '&&' and leftarg = 'geomcollection3D' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_geometry3D, tbl_tgeompoint3D WHERE g @> temp )
WHERE op = '@>' and leftarg = 'geomcollection3D' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_geometry3D, tbl_tgeompoint3D WHERE g <@ temp )
WHERE op = '<@' and leftarg = 'geomcollection3D' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_geometry3D, tbl_tgeompoint3D WHERE g -|- temp )
WHERE op = '-|-' and leftarg = 'geomcollection3D' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_geometry3D, tbl_tgeompoint3D WHERE g ~= temp )
WHERE op = '~=' and leftarg = 'geomcollection3D' and rightarg = 'tgeompoint3D';

UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_timestamptz, tbl_tgeompoint3D WHERE t && temp )
WHERE op = '&&' and leftarg = 'timestamptz' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_timestamptz, tbl_tgeompoint3D WHERE t @> temp )
WHERE op = '@>' and leftarg = 'timestamptz' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_timestamptz, tbl_tgeompoint3D WHERE t <@ temp )
WHERE op = '<@' and leftarg = 'timestamptz' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_timestamptz, tbl_tgeompoint3D WHERE t -|- temp )
WHERE op = '-|-' and leftarg = 'timestamptz' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_timestamptz, tbl_tgeompoint3D WHERE t ~= temp )
WHERE op = '~=' and leftarg = 'timestamptz' and rightarg = 'tgeompoint3D';

UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_timestampset, tbl_tgeompoint3D WHERE ts && temp )
WHERE op = '&&' and leftarg = 'timestampset' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_timestampset, tbl_tgeompoint3D WHERE ts @> temp )
WHERE op = '@>' and leftarg = 'timestampset' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_timestampset, tbl_tgeompoint3D WHERE ts <@ temp )
WHERE op = '<@' and leftarg = 'timestampset' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_timestampset, tbl_tgeompoint3D WHERE ts -|- temp )
WHERE op = '-|-' and leftarg = 'timestampset' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_timestampset, tbl_tgeompoint3D WHERE ts ~= temp )
WHERE op = '~=' and leftarg = 'timestampset' and rightarg = 'tgeompoint3D';

UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_period, tbl_tgeompoint3D WHERE p && temp )
WHERE op = '&&' and leftarg = 'period' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_period, tbl_tgeompoint3D WHERE p @> temp )
WHERE op = '@>' and leftarg = 'period' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_period, tbl_tgeompoint3D WHERE p <@ temp )
WHERE op = '<@' and leftarg = 'period' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_period, tbl_tgeompoint3D WHERE p -|- temp )
WHERE op = '-|-' and leftarg = 'period' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_period, tbl_tgeompoint3D WHERE p ~= temp )
WHERE op = '~=' and leftarg = 'period' and rightarg = 'tgeompoint3D';

UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_periodset, tbl_tgeompoint3D WHERE ps && temp )
WHERE op = '&&' and leftarg = 'periodset' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_periodset, tbl_tgeompoint3D WHERE ps @> temp )
WHERE op = '@>' and leftarg = 'periodset' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_periodset, tbl_tgeompoint3D WHERE ps <@ temp )
WHERE op = '<@' and leftarg = 'periodset' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_periodset, tbl_tgeompoint3D WHERE ps -|- temp )
WHERE op = '-|-' and leftarg = 'periodset' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_periodset, tbl_tgeompoint3D WHERE ps ~= temp )
WHERE op = '~=' and leftarg = 'periodset' and rightarg = 'tgeompoint3D';

UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_stbox, tbl_tgeompoint3D WHERE b && temp )
WHERE op = '&&' and leftarg = 'stbox' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_stbox, tbl_tgeompoint3D WHERE b @> temp )
WHERE op = '@>' and leftarg = 'stbox' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_stbox, tbl_tgeompoint3D WHERE b <@ temp )
WHERE op = '<@' and leftarg = 'stbox' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_stbox, tbl_tgeompoint3D WHERE b -|- temp )
WHERE op = '-|-' and leftarg = 'stbox' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_stbox, tbl_tgeompoint3D WHERE b ~= temp )
WHERE op = '~=' and leftarg = 'stbox' and rightarg = 'tgeompoint3D';

-------------------------------------------------------------------------------
-- <type> op tgeogpoint3D

UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_geography3D, tbl_tgeogpoint3D WHERE g && temp )
WHERE op = '&&' and leftarg = 'geogcollection3D' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_geography3D, tbl_tgeogpoint3D WHERE g @> temp )
WHERE op = '@>' and leftarg = 'geogcollection3D' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_geography3D, tbl_tgeogpoint3D WHERE g <@ temp )
WHERE op = '<@' and leftarg = 'geogcollection3D' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_geography3D, tbl_tgeogpoint3D WHERE g -|- temp )
WHERE op = '-|-' and leftarg = 'geogcollection3D' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_geography3D, tbl_tgeogpoint3D WHERE g ~= temp )
WHERE op = '~=' and leftarg = 'geogcollection3D' and rightarg = 'tgeogpoint3D';

UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_timestamptz, tbl_tgeogpoint3D WHERE t && temp )
WHERE op = '&&' and leftarg = 'timestamptz' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_timestamptz, tbl_tgeogpoint3D WHERE t @> temp )
WHERE op = '@>' and leftarg = 'timestamptz' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_timestamptz, tbl_tgeogpoint3D WHERE t <@ temp )
WHERE op = '<@' and leftarg = 'timestamptz' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_timestamptz, tbl_tgeogpoint3D WHERE t -|- temp )
WHERE op = '-|-' and leftarg = 'timestamptz' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_timestamptz, tbl_tgeogpoint3D WHERE t ~= temp )
WHERE op = '~=' and leftarg = 'timestamptz' and rightarg = 'tgeogpoint3D';

UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_timestampset, tbl_tgeogpoint3D WHERE ts && temp )
WHERE op = '&&' and leftarg = 'timestampset' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_timestampset, tbl_tgeogpoint3D WHERE ts @> temp )
WHERE op = '@>' and leftarg = 'timestampset' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_timestampset, tbl_tgeogpoint3D WHERE ts <@ temp )
WHERE op = '<@' and leftarg = 'timestampset' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_timestampset, tbl_tgeogpoint3D WHERE ts -|- temp )
WHERE op = '-|-' and leftarg = 'timestampset' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_timestampset, tbl_tgeogpoint3D WHERE ts ~= temp )
WHERE op = '~=' and leftarg = 'timestampset' and rightarg = 'tgeogpoint3D';

UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_period, tbl_tgeogpoint3D WHERE p && temp )
WHERE op = '&&' and leftarg = 'period' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_period, tbl_tgeogpoint3D WHERE p @> temp )
WHERE op = '@>' and leftarg = 'period' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_period, tbl_tgeogpoint3D WHERE p <@ temp )
WHERE op = '<@' and leftarg = 'period' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_period, tbl_tgeogpoint3D WHERE p -|- temp )
WHERE op = '-|-' and leftarg = 'period' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_period, tbl_tgeogpoint3D WHERE p ~= temp )
WHERE op = '~=' and leftarg = 'period' and rightarg = 'tgeogpoint3D';

UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_periodset, tbl_tgeogpoint3D WHERE ps && temp )
WHERE op = '&&' and leftarg = 'periodset' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_periodset, tbl_tgeogpoint3D WHERE ps @> temp )
WHERE op = '@>' and leftarg = 'periodset' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_periodset, tbl_tgeogpoint3D WHERE ps <@ temp )
WHERE op = '<@' and leftarg = 'periodset' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_periodset, tbl_tgeogpoint3D WHERE ps -|- temp )
WHERE op = '-|-' and leftarg = 'periodset' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_periodset, tbl_tgeogpoint3D WHERE ps ~= temp )
WHERE op = '~=' and leftarg = 'periodset' and rightarg = 'tgeogpoint3D';

UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_geodstbox, tbl_tgeogpoint3D WHERE b && temp )
WHERE op = '&&' and leftarg = 'stbox' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_geodstbox, tbl_tgeogpoint3D WHERE b @> temp )
WHERE op = '@>' and leftarg = 'stbox' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_geodstbox, tbl_tgeogpoint3D WHERE b <@ temp )
WHERE op = '<@' and leftarg = 'stbox' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_geodstbox, tbl_tgeogpoint3D WHERE b -|- temp )
WHERE op = '-|-' and leftarg = 'stbox' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_geodstbox, tbl_tgeogpoint3D WHERE b ~= temp )
WHERE op = '~=' and leftarg = 'stbox' and rightarg = 'tgeogpoint3D';

-------------------------------------------------------------------------------
-- tgeompoint3D op <type>

UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_geometry3D WHERE temp && g )
WHERE op = '&&' and leftarg = 'tgeompoint3D' and rightarg = 'geomcollection3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_geometry3D WHERE temp @> g )
WHERE op = '@>' and leftarg = 'tgeompoint3D' and rightarg = 'geomcollection3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_geometry3D WHERE temp <@ g )
WHERE op = '<@' and leftarg = 'tgeompoint3D' and rightarg = 'geomcollection3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_geometry3D WHERE temp -|- g )
WHERE op = '-|-' and leftarg = 'tgeompoint3D' and rightarg = 'geomcollection3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_geometry3D WHERE temp ~= g )
WHERE op = '~=' and leftarg = 'tgeompoint3D' and rightarg = 'geomcollection3D';

UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_timestamptz WHERE temp && t )
WHERE op = '&&' and leftarg = 'tgeompoint3D' and rightarg = 'timestamptz';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_timestamptz WHERE temp @> t )
WHERE op = '@>' and leftarg = 'tgeompoint3D' and rightarg = 'timestamptz';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_timestamptz WHERE temp <@ t )
WHERE op = '<@' and leftarg = 'tgeompoint3D' and rightarg = 'timestamptz';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_timestamptz WHERE temp -|- t )
WHERE op = '-|-' and leftarg = 'tgeompoint3D' and rightarg = 'timestamptz';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_timestamptz WHERE temp ~= t )
WHERE op = '~=' and leftarg = 'tgeompoint3D' and rightarg = 'timestamptz';

UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_timestampset WHERE temp && ts )
WHERE op = '&&' and leftarg = 'tgeompoint3D' and rightarg = 'timestampset';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_timestampset WHERE temp @> ts )
WHERE op = '@>' and leftarg = 'tgeompoint3D' and rightarg = 'timestampset';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_timestampset WHERE temp <@ ts )
WHERE op = '<@' and leftarg = 'tgeompoint3D' and rightarg = 'timestampset';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_timestampset WHERE temp -|- ts )
WHERE op = '-|-' and leftarg = 'tgeompoint3D' and rightarg = 'timestampset';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_timestampset WHERE temp ~= ts )
WHERE op = '~=' and leftarg = 'tgeompoint3D' and rightarg = 'timestampset';

UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_period WHERE temp && p )
WHERE op = '&&' and leftarg = 'tgeompoint3D' and rightarg = 'period';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_period WHERE temp @> p )
WHERE op = '@>' and leftarg = 'tgeompoint3D' and rightarg = 'period';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_period WHERE temp <@ p )
WHERE op = '<@' and leftarg = 'tgeompoint3D' and rightarg = 'period';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_period WHERE temp -|- p )
WHERE op = '-|-' and leftarg = 'tgeompoint3D' and rightarg = 'period';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_period WHERE temp ~= p )
WHERE op = '~=' and leftarg = 'tgeompoint3D' and rightarg = 'period';

UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_periodset WHERE temp && ps )
WHERE op = '&&' and leftarg = 'tgeompoint3D' and rightarg = 'periodset';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_periodset WHERE temp @> ps )
WHERE op = '@>' and leftarg = 'tgeompoint3D' and rightarg = 'periodset';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_periodset WHERE temp <@ ps )
WHERE op = '<@' and leftarg = 'tgeompoint3D' and rightarg = 'periodset';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_periodset WHERE temp -|- ps )
WHERE op = '-|-' and leftarg = 'tgeompoint3D' and rightarg = 'periodset';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_periodset WHERE temp ~= ps )
WHERE op = '~=' and leftarg = 'tgeompoint3D' and rightarg = 'periodset';

UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_stbox WHERE temp && b )
WHERE op = '&&' and leftarg = 'tgeompoint3D' and rightarg = 'stbox';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_stbox WHERE temp @> b )
WHERE op = '@>' and leftarg = 'tgeompoint3D' and rightarg = 'stbox';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_stbox WHERE temp <@ b )
WHERE op = '<@' and leftarg = 'tgeompoint3D' and rightarg = 'stbox';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_stbox WHERE temp -|- b )
WHERE op = '-|-' and leftarg = 'tgeompoint3D' and rightarg = 'stbox';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D, tbl_stbox WHERE temp ~= b )
WHERE op = '~=' and leftarg = 'tgeompoint3D' and rightarg = 'stbox';

UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D t1, tbl_tgeompoint3D t2 WHERE t1.temp && t2.temp )
WHERE op = '&&' and leftarg = 'tgeompoint3D' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D t1, tbl_tgeompoint3D t2 WHERE t1.temp @> t2.temp )
WHERE op = '@>' and leftarg = 'tgeompoint3D' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D t1, tbl_tgeompoint3D t2 WHERE t1.temp <@ t2.temp )
WHERE op = '<@' and leftarg = 'tgeompoint3D' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D t1, tbl_tgeompoint3D t2 WHERE t1.temp -|- t2.temp )
WHERE op = '-|-' and leftarg = 'tgeompoint3D' and rightarg = 'tgeompoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeompoint3D t1, tbl_tgeompoint3D t2 WHERE t1.temp ~= t2.temp )
WHERE op = '~=' and leftarg = 'tgeompoint3D' and rightarg = 'tgeompoint3D';

-------------------------------------------------------------------------------
-- tgeogpoint3D op <type>

UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_geography3D WHERE temp && g )
WHERE op = '&&' and leftarg = 'tgeogpoint3D' and rightarg = 'geogcollection3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_geography3D WHERE temp @> g )
WHERE op = '@>' and leftarg = 'tgeogpoint3D' and rightarg = 'geogcollection3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_geography3D WHERE temp <@ g )
WHERE op = '<@' and leftarg = 'tgeogpoint3D' and rightarg = 'geogcollection3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_geography3D WHERE temp -|- g )
WHERE op = '-|-' and leftarg = 'tgeogpoint3D' and rightarg = 'geogcollection3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_geography3D WHERE temp ~= g )
WHERE op = '~=' and leftarg = 'tgeogpoint3D' and rightarg = 'geogcollection3D';

UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_timestamptz WHERE temp && t )
WHERE op = '&&' and leftarg = 'tgeogpoint3D' and rightarg = 'timestamptz';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_timestamptz WHERE temp @> t )
WHERE op = '@>' and leftarg = 'tgeogpoint3D' and rightarg = 'timestamptz';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_timestamptz WHERE temp <@ t )
WHERE op = '<@' and leftarg = 'tgeogpoint3D' and rightarg = 'timestamptz';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_timestamptz WHERE temp -|- t )
WHERE op = '-|-' and leftarg = 'tgeogpoint3D' and rightarg = 'timestamptz';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_timestamptz WHERE temp ~= t )
WHERE op = '~=' and leftarg = 'tgeogpoint3D' and rightarg = 'timestamptz';

UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_timestampset WHERE temp && ts )
WHERE op = '&&' and leftarg = 'tgeogpoint3D' and rightarg = 'timestampset';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_timestampset WHERE temp @> ts )
WHERE op = '@>' and leftarg = 'tgeogpoint3D' and rightarg = 'timestampset';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_timestampset WHERE temp <@ ts )
WHERE op = '<@' and leftarg = 'tgeogpoint3D' and rightarg = 'timestampset';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_timestampset WHERE temp -|- ts )
WHERE op = '-|-' and leftarg = 'tgeogpoint3D' and rightarg = 'timestampset';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_timestampset WHERE temp ~= ts )
WHERE op = '~=' and leftarg = 'tgeogpoint3D' and rightarg = 'timestampset';

UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_period WHERE temp && p )
WHERE op = '&&' and leftarg = 'tgeogpoint3D' and rightarg = 'period';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_period WHERE temp @> p )
WHERE op = '@>' and leftarg = 'tgeogpoint3D' and rightarg = 'period';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_period WHERE temp <@ p )
WHERE op = '<@' and leftarg = 'tgeogpoint3D' and rightarg = 'period';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_period WHERE temp -|- p )
WHERE op = '-|-' and leftarg = 'tgeogpoint3D' and rightarg = 'period';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_period WHERE temp ~= p )
WHERE op = '~=' and leftarg = 'tgeogpoint3D' and rightarg = 'period';

UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_periodset WHERE temp && ps )
WHERE op = '&&' and leftarg = 'tgeogpoint3D' and rightarg = 'periodset';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_periodset WHERE temp @> ps )
WHERE op = '@>' and leftarg = 'tgeogpoint3D' and rightarg = 'periodset';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_periodset WHERE temp <@ ps )
WHERE op = '<@' and leftarg = 'tgeogpoint3D' and rightarg = 'periodset';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_periodset WHERE temp -|- ps )
WHERE op = '-|-' and leftarg = 'tgeogpoint3D' and rightarg = 'periodset';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_periodset WHERE temp ~= ps )
WHERE op = '~=' and leftarg = 'tgeogpoint3D' and rightarg = 'periodset';

UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_geodstbox WHERE temp && b )
WHERE op = '&&' and leftarg = 'tgeogpoint3D' and rightarg = 'stbox';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_geodstbox WHERE temp @> b )
WHERE op = '@>' and leftarg = 'tgeogpoint3D' and rightarg = 'stbox';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_geodstbox WHERE temp <@ b )
WHERE op = '<@' and leftarg = 'tgeogpoint3D' and rightarg = 'stbox';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_geodstbox WHERE temp -|- b )
WHERE op = '-|-' and leftarg = 'tgeogpoint3D' and rightarg = 'stbox';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D, tbl_geodstbox WHERE temp ~= b )
WHERE op = '~=' and leftarg = 'tgeogpoint3D' and rightarg = 'stbox';

UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D t1, tbl_tgeogpoint3D t2 WHERE t1.temp && t2.temp )
WHERE op = '&&' and leftarg = 'tgeogpoint3D' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D t1, tbl_tgeogpoint3D t2 WHERE t1.temp @> t2.temp )
WHERE op = '@>' and leftarg = 'tgeogpoint3D' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D t1, tbl_tgeogpoint3D t2 WHERE t1.temp <@ t2.temp )
WHERE op = '<@' and leftarg = 'tgeogpoint3D' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D t1, tbl_tgeogpoint3D t2 WHERE t1.temp -|- t2.temp )
WHERE op = '-|-' and leftarg = 'tgeogpoint3D' and rightarg = 'tgeogpoint3D';
UPDATE test_geoboundboxops3d
SET spgistidx = ( SELECT count(*) FROM tbl_tgeogpoint3D t1, tbl_tgeogpoint3D t2 WHERE t1.temp ~= t2.temp )
WHERE op = '~=' and leftarg = 'tgeogpoint3D' and rightarg = 'tgeogpoint3D';

-------------------------------------------------------------------------------

DROP INDEX IF EXISTS tbl_tgeompoint3D_spgist_idx;
DROP INDEX IF EXISTS tbl_tgeogpoint3D_spgist_idx;

-------------------------------------------------------------------------------

SELECT * FROM test_geoboundboxops3d
WHERE noidx <> spgistidx
ORDER BY op, leftarg, rightarg;

DROP TABLE test_geoboundboxops3d;

-------------------------------------------------------------------------------
