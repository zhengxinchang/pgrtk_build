impl :: bincode :: Decode for Fragment
{
    fn decode < __D : :: bincode :: de :: Decoder > (decoder : & mut __D)
    ->core :: result :: Result < Self, :: bincode :: error :: DecodeError >
    {
        let variant_index = < u32 as :: bincode :: Decode >:: decode(decoder)
        ?; match variant_index
        {
            0u32
            =>Ok(Self ::AlnSegments
            { 0 : :: bincode :: Decode :: decode(decoder) ?, }), 1u32
            =>Ok(Self ::Prefix
            { 0 : :: bincode :: Decode :: decode(decoder) ?, }), 2u32
            =>Ok(Self ::Internal
            { 0 : :: bincode :: Decode :: decode(decoder) ?, }), 3u32
            =>Ok(Self ::Suffix
            { 0 : :: bincode :: Decode :: decode(decoder) ?, }), variant
            =>Err(:: bincode :: error :: DecodeError :: UnexpectedVariant
            {
                found : variant, type_name : "Fragment", allowed : &:: bincode
                :: error :: AllowedEnumVariants :: Range { min: 0, max: 3 }
            })
        }
    }
} impl < '__de > :: bincode :: BorrowDecode < '__de > for Fragment
{
    fn borrow_decode < __D : :: bincode :: de :: BorrowDecoder < '__de > >
    (decoder : & mut __D) ->core :: result :: Result < Self, :: bincode ::
    error :: DecodeError >
    {
        let variant_index = < u32 as :: bincode :: Decode >:: decode(decoder)
        ?; match variant_index
        {
            0u32
            =>Ok(Self ::AlnSegments
            { 0 : :: bincode :: BorrowDecode :: borrow_decode(decoder) ?, }),
            1u32
            =>Ok(Self ::Prefix
            { 0 : :: bincode :: BorrowDecode :: borrow_decode(decoder) ?, }),
            2u32
            =>Ok(Self ::Internal
            { 0 : :: bincode :: BorrowDecode :: borrow_decode(decoder) ?, }),
            3u32
            =>Ok(Self ::Suffix
            { 0 : :: bincode :: BorrowDecode :: borrow_decode(decoder) ?, }),
            variant
            =>Err(:: bincode :: error :: DecodeError :: UnexpectedVariant
            {
                found : variant, type_name : "Fragment", allowed : &:: bincode
                :: error :: AllowedEnumVariants :: Range { min: 0, max: 3 }
            })
        }
    }
}