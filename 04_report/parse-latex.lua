if FORMAT == 'latex' then
  return {}
end

function RawBlock (raw)
  if raw.format:match 'tex' then
    print(raw.text)
    return pandoc.read(raw.text, 'latex').blocks
  end
end