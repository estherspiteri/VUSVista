import { IPublicationPreview } from "../../models/publication-view.model";
import { IVUSSummary } from "../../models/vus-summary.model";

export interface IGetPublicationsByVariantIdRequest {
  variantId: string;
}

export interface IGetPublicationsByVariantIdResponse {
  isSuccess: boolean;
  publications?: IPublicationPreview[];
  variant: IVUSSummary;
}

export interface IGetPublicationsByRsidAndWithOptionalTextRequest {
  variantId: string;
  rsid?: string | null;
  optionalText: string;
}

export interface IGetPublicationsByRsidAndWithOptionalTextResponse {
  isSuccess: boolean;
  publications?: IPublicationPreview[];
  variant: IVUSSummary;
}
