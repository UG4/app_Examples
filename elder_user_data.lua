--------------------------------
-- User Data Functions 2D
--------------------------------

if dim == 2 then 

function ConcentrationStart(x, y, t)
  if y == 150 then
    if x > 150 and x < 450 then
    return 1.0
    end
  end
  return 0.0
end

function PressureStart(x, y, t)
  return 9810 * (150 - y)
end

function ConcentrationDirichletBnd(x, y, t)
  if y == 150 then
    if x > 150 and x < 450 then
      return true, 1.0
    end
  end
  if y == 0.0 then
    return true, 0.0
  end

  return false, 0.0
end

function PressureDirichletBnd(x, y, t)
  if y == 150 then
    if x == 0.0 or x == 600 then
      return true, 9810 * (150 - y)
    end
  end
  
  return false, 0.0
end

function Porosity(x,y,t)
  return 0.1
end

else 
--------------------------------
-- User Data Functions 3D
--------------------------------
  
function ConcentrationStart(x, y, z, t)
  if z == 150 then
    if y > 150 and y < 450 then
      if x > 150 and x < 450 then
        return 1.0
      end
    end
  end
  return 0.0
end

function PressureStart(x, y, z, t)
  return 9810 * (150 - z)
end

function ConcentrationDirichletBnd(x, y, z, t)
  if z == 150 then
    if y > 150 and y < 450 then
      if x > 150 and x < 450 then
        return true, 1.0
      end
    end
  end
  if z == 0.0 then
    return true, 0.0
  end
  
  return false, 0.0
end

function PressureDirichletBnd(x, y, z, t)
  if z == 150 then
    if y == 0.0 or y == 600 then
      if x == 0.0 or x == 600 then
        return true, 9810 * (150 - z)
      end
    end
  end
  
  return false, 0.0
end

function Porosity(x,y,z,t)
  return 0.1
end

end
