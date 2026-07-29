-- eq-number-format.lua
-- Changes the inline equation label Quarto embeds in display math from
--   \qquad(N)  →  \qquad(\text{Eq. }N)
--
-- Background: for non-HTML/non-LaTeX output (i.e. docx) Quarto's crossref
-- filter appends  \qquad(N)  to the raw LaTeX of every labelled display-math
-- block (see eqQquad() in crossref.lua).  Pandoc then converts the whole
-- expression — equation + number — to OMML, so the "(N)" appears as a
-- mathematical expression inside the equation box, not as margin text.
-- The \qquad provides separation, but (N) still reads as maths.
--
-- This filter (which runs *after* the crossref filter) replaces the pattern
-- so the label renders as upright text, clearly marking it as a reference
-- number rather than a mathematical subexpression.
--
-- Only applied to docx output; other formats use different code paths.

if FORMAT ~= "docx" then
  return {}
end

function Math(el)
  if el.mathtype ~= "DisplayMath" then
    return nil
  end
  -- Match \qquad(…) at the very end of the LaTeX string.
  -- The number may be a bare digit ("1") or section-qualified ("2.1").
  local new_text, n = el.text:gsub(
    "\\qquad%(([^%)]+)%)%s*$",
    "\\qquad(\\text{Eq.\\ }%1)"
  )
  if n > 0 then
    el.text = new_text
    return el
  end
end
