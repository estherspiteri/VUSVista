import { IPublicationPreview, IVUSSummary } from "../../models/publication-view.model";

export interface IGetPublicationsByVariantIdRequest {
  variantId: string;
}

export interface IGetPublicationsByVariantIdResponse {
  isSuccess: boolean;
  publications?: IPublicationPreview[];
  variant: IVUSSummary;
}
