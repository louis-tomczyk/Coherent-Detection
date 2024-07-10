function dlambda = dnu2dlambda(dnu,nu)

  % ---------------------------------------------
  % ----- INFORMATIONS -----
  %   Function name   : 
  %   Author          : louis tomczyk
  %   Institution     : Telecom Paris
  %   Email           : louis.tomczyk.work@gmail.com
  %   Date            : 2022-08-23
  %   Version         : 1.0
  %
  % ----- MAIN IDEA -----
  % ----- INPUTS -----
  %   DNU [GHz]   linewidth, typically - 37.5 for channel spacing (C-band)
  %                                    - 12.5 for OSA resolution 
  %                         both @1550 [nm]
  %   NU  [THz]   frequency, typically - 193.41 [THz] == 1550 [nm] 
  %                                    - 230.61 [THz] == 1300 [nm]
  %                                    - 281.76 [THz] == 1064 [nm]
  %                                    - 473.76 [THz] == 632.8[nm]
  %                                    - 508.47 [THz] == 589.6[nm]
  %                                    - 563.52 [THz] == 532  [nm]
  %
  % ----- BIBLIOGRAPHY -----
  % ---------------------------------------------
  
  
    dnu     = dnu*1e9;            % [Hz]
    nu      = nu*1e12;            % [Hz]
  
    c       = 299792458;
    dlambda = c.*dnu/(nu.^2)*1e9; % [nm]
  
  end
  
  