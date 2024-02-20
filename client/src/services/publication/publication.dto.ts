import {
  IPublicationPreview,
} from "../../models/publication-search.model";

export interface IGetPublicationsByVariantIdRequest {
  variantId: string;
}

export interface IGetPublicationsByVariantIdResponse {
  isSuccess: boolean;
  publications?: IPublicationPreview[];
}
